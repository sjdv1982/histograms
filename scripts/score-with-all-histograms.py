# scores with all atomtype-atomtype interactions
# Assumes that <complex-dir>/receptor-reduced exists
# as well as <complex-dir>/coordinates/X-Y.npy, where X and Y are receptor and ligand atomtypes
# histogram_dir must contain X-Y.json files.

import sys, os
import numpy as np
import json

complex_dir = sys.argv[1]
motif = sys.argv[2]
histogram_dir = sys.argv[3]

receptor_pdb = os.path.join(complex_dir, "receptor-reduced.pdb")
assert os.path.exists(receptor_pdb)

receptor = {}
for l in open(complex_dir + "/receptor-reduced.pdb"):
    if not l.startswith("ATOM"):
        continue
    atomtype = int(l[57:59])
    x, y, z = float(l[30:38]),float(l[38:46]),float(l[46:54])
    if atomtype not in receptor:
        receptor[atomtype] = []
    receptor[atomtype].append((x,y,z))
for atomtype in receptor:
    receptor[atomtype] = np.array(receptor[atomtype])

def score_with_histogram(histogram, receptor_atomtype, ligand_coordinates):
    nstruc = len(ligand_coordinates)
    scores = np.zeros(nstruc)
    receptor_coordinates = receptor[receptor_atomtype]

    for coor in receptor_coordinates:
        dist = np.linalg.norm(ligand_coordinates - coor, axis=2)
        for n in range(len(histogram["rank_chunks"])):
            low_rank, high_rank = histogram["rank_chunks"][n]
            distance_bins = histogram["distance_bins"][n]
            distance_thresholds = np.array([0] + [v[0] for v in distance_bins])
            rank_scores = np.array([v[1] for v in distance_bins] + [0])
            for lig_coor in range(ligand_coordinates.shape[1]):
                curr_dist = dist[low_rank:high_rank, lig_coor]
                curr_bin_dist = np.digitize(curr_dist, distance_thresholds) - 1
                curr_scores = rank_scores[curr_bin_dist]
                sc = scores[low_rank:high_rank]
                sc += curr_scores
    return scores

nstruc = None

all_ligand_coordinates = {}
for ligand_atomtype in range(1, 99):
    f = "{}-{}.npy".format(motif, ligand_atomtype)
    ligand_coordinates_file = os.path.join(complex_dir, "coordinates", f)
    if not os.path.exists(ligand_coordinates_file):
        continue
    ligand_coordinates = np.load(ligand_coordinates_file)
    if nstruc is None:
        nstruc = len(ligand_coordinates)        
    else:
        assert len(ligand_coordinates) == nstruc, ligand_coordinates_file        
    all_ligand_coordinates[ligand_atomtype] = ligand_coordinates

jobs = []
for ligand_atomtype in range(1, 99):
    if ligand_atomtype not in all_ligand_coordinates:
        continue
    ligand_coordinates = all_ligand_coordinates[ligand_atomtype]
    for receptor_atomtype in range(1, 99):
        f = "{}-{}.json".format(receptor_atomtype, ligand_atomtype)
        histogram_file = os.path.join(histogram_dir, f)
        if not os.path.exists(histogram_file):
            continue
        histogram = json.load(open(histogram_file))
        assert len(histogram["rank_chunks"]) == len(histogram["distance_bins"])        
        jobs.append((histogram, receptor_atomtype, ligand_coordinates))    

"""
# non-parallel
scores = np.zeros(nstruc)
for job in jobs:
    histogram, receptor_atomtype, ligand_coordinates = job
    scores += score_with_histogram(histogram, receptor_atomtype, ligand_coordinates)
"""

# parallel
import multiprocessing
pool = multiprocessing.Pool()
all_scores = pool.starmap(score_with_histogram, jobs)
scores = np.sum(all_scores, axis=0)

for score in scores:
    print("%.6f" % score)
# scores with a single atomtype-atomtype interaction

import sys, os
import numpy as np
import json

receptor_file = sys.argv[1] # e.g. 1A9N/receptor-reduced.pdb
receptor_atomtype = int(sys.argv[2])  # e.g. 6
histogram_file = sys.argv[3] # e.g. histograms/6-38.json
ligand_coordinates_file = sys.argv[4]  #e.g 1A9N/coordinates/AUU-38.npy

histogram = json.load(open(histogram_file))
assert len(histogram["rank_chunks"]) == len(histogram["distance_bins"])

receptor_coordinates = []
for l in open(receptor_file):
    if not l.startswith("ATOM"):
        continue
    atomtype = int(l[57:59])
    if atomtype != receptor_atomtype:
        continue
    x, y, z = float(l[30:38]),float(l[38:46]),float(l[46:54])
    receptor_coordinates.append((x,y,z))
receptor_coordinates = np.array(receptor_coordinates)

ligand_coordinates = np.load(ligand_coordinates_file)
nstruc = len(ligand_coordinates)
scores = np.zeros(nstruc)
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

for score in scores:
    print("%.6f" % score)
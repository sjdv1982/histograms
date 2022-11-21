import sys, os
import numpy as np
import random

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from score_with_histograms_discrete import calc_score, load_histograms, load_receptor_coordinates

receptor_file = sys.argv[1] # e.g. 1A9N/receptor-reduced.pdb
coordinates_dir = sys.argv[2] # must contain X-Y.npy, where X is the motif and Y is the ligand atom type
motif = sys.argv[3]  # e.g. GUG
histogram_dir = sys.argv[4] # must contain Z-Y.json, where Z is the receptor atom type and Y is the ligand atom type
cache_dir = sys.argv[5]  # must contain discrete-Y-Q.npy (generated with discretize-coordinates.py), where Y is the ligand atom type and Q is the checksum of the rank chunk scheme

all_receptor_coordinates = load_receptor_coordinates(receptor_file)
ligand_atoms = []
for ligand_atomtype in range(32, 48+1):
    ligand_file = os.path.join(coordinates_dir, "{}-{}.npy".format(motif, ligand_atomtype))
    if os.path.exists(ligand_file):
        ligand_atoms.append((ligand_atomtype, ligand_file))

random.shuffle(ligand_atoms)

receptor_atomtypes = sorted(list(all_receptor_coordinates.keys()))
scores = None
for ligand_atomtype, ligand_file in ligand_atoms:
    histogram_files = ["{}/{}-{}.json".format(histogram_dir, k, ligand_atomtype) for k in receptor_atomtypes]
    ligand_coordinates = np.load(ligand_file)
    histograms = load_histograms(histogram_files, ligand_atomtype, len(ligand_coordinates))
    curr_scores = calc_score(all_receptor_coordinates, ligand_coordinates, ligand_atomtype, histograms, cache_dir)
    if scores is None:
        scores = curr_scores
    else:
        scores += curr_scores

for score in scores:
    print("%.6f" % score)

import sys, os
import numpy as np
import random
import glob

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from score_with_histograms_discrete import calc_score, load_histograms, load_receptor_coordinates

receptor_file = sys.argv[1] # e.g. 1A9N/receptor-reduced.pdb
motif = sys.argv[2]  # e.g. GUG
nstruc = int(sys.argv[3]) # the number of ligand structures
histogram_dir = sys.argv[4] # must contain Z-Y.json, where Z is the receptor atom type and Y is the ligand atom type
cache_dir = sys.argv[5]  # must contain discrete-Y-Q.npy (generated with discretize-coordinates.py), where Y is the ligand atom type and Q is the checksum of the rank chunk scheme
                         # The cache dir contents determine what ligand atomtypes are actually there

all_receptor_coordinates = load_receptor_coordinates(receptor_file)
ligand_atomtypes = []
for ligand_atomtype in range(32, 48+1):
    if len(glob.glob("{}/discrete-{}-*".format(cache_dir, ligand_atomtype))):
        ligand_atomtypes.append(ligand_atomtype)

if not ligand_atomtypes:
    raise Exception("Nothing found in the cache dir")

random.shuffle(ligand_atomtypes)

receptor_atomtypes = sorted(list(all_receptor_coordinates.keys()))
scores = None
for ligand_atomtype in ligand_atomtypes:
    histogram_files0 = ["{}/{}-{}.json".format(histogram_dir, k, ligand_atomtype) for k in receptor_atomtypes]
    histogram_files = []
    for f in histogram_files0:
        if not os.path.exists(f):
            print("WARNING: {} does not exist".format(f),file=sys.stderr)
        else:
            histogram_files.append(f)
    histograms = load_histograms(histogram_files, ligand_atomtype, nstruc)
    curr_scores = calc_score(all_receptor_coordinates, None, ligand_atomtype, histograms, cache_dir, nstruc=nstruc)
    if scores is None:
        scores = curr_scores
    else:
        scores += curr_scores

for score in scores:
    print("%.6f" % score)

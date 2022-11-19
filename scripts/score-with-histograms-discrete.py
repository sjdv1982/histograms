# scores with a single atomtype-atomtype interaction
# use discretization at 1A with linear interpolation, if it makes it faster

import sys, os
import numpy as np
import json
from hashlib import sha3_256
import pandas as pd

def discretize_coordinates(coordinates, rank_chunks):
    bits10 = np.uint32(2**10-1)

    coor_min = np.floor(coordinates.min(axis=0).min(axis=0))
    coordinates = coordinates - coor_min
    rank_start = 0

    result_dtype = np.dtype([("index", np.uint32),("weight", np.uint8)],align=False)
    digit_indices = np.empty(coordinates.shape[:2] + (8,), dtype=result_dtype)

    all_digit_coors = []
    for rank_start, rank_end in rank_chunks:
        chunk = coordinates[rank_start:rank_end]    
        if not len(chunk):
            continue
        discrete = []
        for xyz in range(3):
            coor = chunk[:, :, xyz]
            minus = np.floor(coor).astype(np.uint32)
            plus = minus + 1
            assert plus.max() < 1024
            wminus = plus - coor
            wplus = coor - minus
            discrete.append(((minus, wminus),(plus, wplus)))
        # digit_chunk = np.empty(chunk.shape, np.uint32)
        digit_hash = np.empty((8,) + chunk.shape[:2], np.uint32)
        weights = np.empty((8,) + chunk.shape[:2])
        ind_xyz = 0
        for dcoorx, weightx in discrete[0]:
            # digit_chunk[:, :, 0] = dcoorx
            hash_x = (dcoorx << 20)
            for dcoory, weighty in discrete[1]:
                # digit_chunk[:, :, 1] = dcoory
                hash_y = (dcoory << 10)
                weightxy = weightx * weighty
                for dcoorz, weightz in discrete[2]:
                    # digit_chunk[:, :, 2] = dcoorz
                    digit_hash[ind_xyz] = hash_x + hash_y + dcoorz
                    weights[ind_xyz] = weightxy * weightz
                    ind_xyz += 1
        # slow:
        #uniq, digit_index = np.unique(digit_hash, return_inverse=True) 
        # faster:
        uniq = np.sort(pd.unique(digit_hash.flatten()))
        digit_index = np.searchsorted(uniq, digit_hash) # two-thirds of the computation...

        digit_coors = np.empty((len(uniq), 3))
        digit_coors[:, 0] = (uniq >> 20) & bits10
        digit_coors[:, 1] = (uniq >> 10) & bits10
        digit_coors[:, 2] = uniq & bits10
        digit_index = digit_index.reshape(digit_hash.shape).astype(np.uint32)
        print("Number of atoms: {}, unique on grid: {}".format(chunk.shape[0] * chunk.shape[1], len(uniq)), file=sys.stderr)
        digit_indices["index"][rank_start:rank_end] = np.moveaxis(digit_index, 0, -1)
        weights = np.round(weights * 255).astype(np.uint8)
        digit_indices["weight"][rank_start:rank_end] = np.moveaxis(weights, 0, -1)
        all_digit_coors.append(digit_coors)
    
    digit_coors = np.concatenate(all_digit_coors) + coor_min
    return digit_coors, digit_indices

receptor_file = sys.argv[1] # e.g. 1A9N/receptor-reduced.pdb

ligand_coordinates = sys.argv[2] # e.g GUG-32.npy

ligand_atomtype = int(sys.argv[3]) # e.g. 32

histogram_files = sys.argv[4:] # e.g. histograms/6-38.json. 
    # Can be a list of Y-X.json, with the same X (ligand atom type)

for h in histogram_files:
    hh = os.path.splitext(os.path.split(h)[1])[0]
    _, lig0 = [int(k) for k in hh.split("-")]
    if lig0 != ligand_atomtype:
        raise Exception("Incorrect ligand type in histogram files")

all_receptor_coordinates = {}
for l in open(receptor_file):
    if not l.startswith("ATOM"):
        continue
    atomtype = int(l[57:59])
    rec_coor = all_receptor_coordinates.get(atomtype)
    if rec_coor is None:
        rec_coor = []
        all_receptor_coordinates[atomtype] = rec_coor
    x, y, z = float(l[30:38]),float(l[38:46]),float(l[46:54])
    rec_coor.append((x,y,z))


ligand_coordinates = np.load(ligand_coordinates)
nstruc = len(ligand_coordinates)

scores = np.zeros(nstruc)

histograms = {}

for histogram_file in histogram_files:
    with open(histogram_file) as f:
        histogram = json.load(f)
    assert len(histogram["rank_chunks"]) == len(histogram["distance_bins"])

    rank_chunks0 = histogram["rank_chunks"]
    rank_chunks = []
    for rank_start, rank_end in rank_chunks0:
        rank_start, rank_end = int(rank_start), int(rank_end)
        if rank_start >= nstruc:
            continue
        if rank_end > nstruc:
            rank_end = nstruc
        rank_chunks.append((rank_start, rank_end))
    rankchunk_checksum = sha3_256(json.dumps(rank_chunks, indent=2).encode()).digest().hex()
    histograms[histogram_file] = rankchunk_checksum, rank_chunks, histogram

rankchunk_checksum = None
for histogram_file in sorted(histograms, key=lambda k: histograms[k][0]):
    print(histogram_file, file=sys.stderr)
    hh = os.path.splitext(os.path.split(histogram_file)[1])[0]
    receptor_atomtype, _ = [int(k) for k in hh.split("-")]
    if receptor_atomtype not in all_receptor_coordinates:
        continue
    receptor_coordinates = np.array(all_receptor_coordinates[receptor_atomtype])

    rankchunk_checksum0, rank_chunks, histogram = histograms[histogram_file]
    
    if rankchunk_checksum0 != rankchunk_checksum:
        rankchunk_checksum = rankchunk_checksum0

        digit_coordinates, digit_indices = discretize_coordinates(ligand_coordinates, rank_chunks)

    index_pos = 0
    for n, (low_rank, high_rank) in enumerate(rank_chunks):
        curr_indices = digit_indices["index"][low_rank:high_rank]
        curr_weights = digit_indices["weight"][low_rank:high_rank]
        ncoors = curr_indices.max() + 1
        chunk = digit_coordinates[index_pos:index_pos+ncoors]
        index_pos += ncoors

        distance_bins = histogram["distance_bins"][n]
        distance_thresholds = np.array([0] + [v[0] for v in distance_bins])
        rank_potential = np.array([v[1] for v in distance_bins] + [0])

        chunk_scores = np.zeros(len(chunk))
        for coor in receptor_coordinates:
            dist = np.linalg.norm(chunk - coor, axis=1)
            bin_dist = np.digitize(dist, distance_thresholds) - 1
            curr_scores = rank_potential[bin_dist]
            chunk_scores += curr_scores
        curr_scores_vec = chunk_scores[curr_indices]
        curr_scores = (curr_scores_vec/255 * curr_weights).reshape(len(curr_indices), -1).sum(axis=1)
        scores[low_rank:high_rank] += curr_scores

for score in scores:
    print("%.6f" % score)

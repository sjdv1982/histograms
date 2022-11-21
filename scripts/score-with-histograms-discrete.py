# scores with a single atomtype-atomtype interaction
# use discretization at 1A with linear interpolation, if it makes it faster

import sys, os
import numpy as np
import json
from hashlib import sha3_256

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from discretize_coordinates import discretize_coordinates

def load_receptor_coordinates(receptor_file):
    all_receptor_coordinates = {}
    for l in open(receptor_file):
        if not l.startswith("ATOM"):
            continue
        atomtype = int(l[57:59])
        if atomtype == 99:
            continue
        rec_coor = all_receptor_coordinates.get(atomtype)
        if rec_coor is None:
            rec_coor = []
            all_receptor_coordinates[atomtype] = rec_coor
        x, y, z = float(l[30:38]),float(l[38:46]),float(l[46:54])
        rec_coor.append((x,y,z))
    for k in list(all_receptor_coordinates.keys()):
        all_receptor_coordinates[k] = np.array(all_receptor_coordinates[k])
    return all_receptor_coordinates

def load_histograms(histogram_files, ligand_atomtype, nstruc):
    for h in histogram_files:
        hh = os.path.splitext(os.path.split(h)[1])[0]
        _, lig0 = [int(k) for k in hh.split("-")]
        if lig0 != ligand_atomtype:
            raise Exception("Incorrect ligand type in histogram files")

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
    return histograms

def calc_score(all_receptor_coordinates, ligand_coordinates, ligand_atomtype, histograms, cache_dir, *, nstruc=None):

    if nstruc is None:
        nstruc = len(ligand_coordinates)
    scores = np.zeros(nstruc)


    rankchunk_checksum = None
    for histogram_file in sorted(histograms, key=lambda k: histograms[k][0]):
        print(histogram_file, file=sys.stderr)
        hh = os.path.splitext(os.path.split(histogram_file)[1])[0]
        receptor_atomtype, _ = [int(k) for k in hh.split("-")]
        if receptor_atomtype not in all_receptor_coordinates:
            continue
        receptor_coordinates = all_receptor_coordinates[receptor_atomtype]

        rankchunk_checksum0, rank_chunks, histogram = histograms[histogram_file]
        
        if rankchunk_checksum0 != rankchunk_checksum:
            rankchunk_checksum = rankchunk_checksum0

            digit_coordinates, digit_indices = None, None
            if cache_dir is not None:
                digit_pattern = os.path.join(cache_dir, "discrete-{}-{}".format(ligand_atomtype, rankchunk_checksum))
                digitfile_1 = digit_pattern + "-coordinates.npy"
                digitfile_2 = digit_pattern + "-indices.npy"
                if os.path.exists(digitfile_1):
                    digit_coordinates = np.load(digitfile_1)
                if os.path.exists(digitfile_2):
                    digit_indices = np.load(digitfile_2)

            if digit_coordinates is None or digit_indices is None:
                if ligand_coordinates is None:
                    raise Exception("Cache miss, and no ligand coordinates provided: {}".format(digit_pattern))
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
    return scores

if __name__ == "__main__":
    receptor_file = sys.argv[1] # e.g. 1A9N/receptor-reduced.pdb

    ligand_coordinates = sys.argv[2] # e.g GUG-32.npy

    ligand_atomtype = int(sys.argv[3]) # e.g. 32

    histogram_files = sys.argv[4:] # e.g. histograms/6-38.json. 
        # Can be a list of Y-X.json, with the same X (ligand atom type)

    cache_dir = None
    try:
        pos = histogram_files.index("--cache") 
        cache_dir = histogram_files[pos+1]
        histogram_files = histogram_files[:pos] + histogram_files[pos+2:]
        del pos
    except IndexError:
        pass

    all_receptor_coordinates = load_receptor_coordinates(receptor_file)
    ligand_coordinates = np.load(ligand_coordinates)
    histograms = load_histograms(histogram_files, len(ligand_coordinates))

    scores = calc_score(all_receptor_coordinates, ligand_coordinates, ligand_atomtype, histograms, cache_dir)
    for score in scores:
        print("%.6f" % score)

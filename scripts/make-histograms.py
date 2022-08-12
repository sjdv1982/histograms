"""Makes histograms from distance arrays.
Input: See build-distance-arrays.py for details
One distance array must be present for every atomtype-atomtype combination
Each distance array contains small bins with (nearnative, nonnative) counts.
This is done separately for each chunk of ATTRACT ranks.

make-histograms first joins the rank chunks into bigger chunks. 
In principle, all possible rank chunk joinings are tried
(but with two speed-ups, see below). 
For every rank chunk joining, a step potential consisting of distance bins (with a log-odds score for each bin)
is then automatically derived, by merging small bins until at least a minimum of near-natives are present in
the merged bin.
Therefore, there is a one-to-one relationship between step potential and joining
The best step potential is selected by discriminative power.
For now, that power is computed as 

    abs(log_odds) * nat, where:
        nat is the number of near-natives in the bin
        log_odds is the log-odds score for that bin
    summed over all bins

In terms of speed up:
- The best joining of chunk 1..N is computed as the best joining of 1..K, plus a single chunk K..N,
  considering all values K between 1 and N. 
  The best joining of 1..K will have been computed before, as N is incremented from 1 to len(rank_chunks)
  This is just an optimization (no approximation)
- Very big chunks K..N are not considered, with "very big" defined as at least twice the biggest chunk in the
  best joining encountered so far (the best joining for 1..2, 1..3, ..., 1..N-1)

make-histograms usually runs within minutes, unless:
- the number of rank chunks is very high (=> many possible joinings)
- minimum_near_natives is very low (=> many bins)

"""
import sys, os
import numpy as np
from math import log
import argparse
import json

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("distance_array_dir",
    help="Name of the directory containing distance arrays, with names such as 1-32.npy"
)
p.add_argument("minimum_near_natives", type=int,
    help="The minimum number of near-native structures in each histogram bin"
)

p.add_argument("output_dir",
    help="The directory where to store the generated histograms, with names such as 1-32.json "
)

p.add_argument("--nparallel", type=int, default=None, help="Number of arrays to compute in parallel. By default, all cores are used.")

args = p.parse_args()

minimum_near_natives = args.minimum_near_natives

def join(disc_dists):
    result = disc_dists[0].copy()
    result[0] = disc_dists[:, 0].sum()
    result[1] = disc_dists[:, 1].sum()
    return result

def calc_score(disc_dist, distance_ranges):
    pos = 0
    score = 0
    near_natives = disc_dist[:, 0]
    non_natives = disc_dist[:, 1]
    tot_near_natives = near_natives.sum()
    if not tot_near_natives:
        return (), 0
    baseline = log(tot_near_natives) - log(non_natives.sum())
    bins, potential = [], []
    while 1:
        nat, nonnat = 0, 0
        while 1:
            if pos == len(disc_dist):
                return list(zip(bins, potential)), score
            nat += near_natives[pos]
            nonnat += non_natives[pos]
            pos += 1
            if nat >= minimum_near_natives:
                break
        log_odds = log(nat) - log(nonnat) - baseline
        potential.append(log_odds)
        score += abs(log_odds) * nat
        if pos == len(distance_ranges):
            bins.append(999)
        else:
            bins.append(float(distance_ranges[pos]))

def make_histogram(dist_array_file, outfile):
    print(dist_array_file)
    dist_array = np.load(dist_array_file) # distance array, generated with build-distance-arrays.py
    distance_ranges = [0] + dist_array["distance_ranges"].tolist()
    rank_ranges = [0] + dist_array["rank_ranges"].tolist()
    data = dist_array["data"]  # rank chunks
    assert data.shape[0] == len(rank_ranges)-1, (data.shape, len(rank_ranges)-1)
    assert data.shape[1] == len(distance_ranges), (data.shape, len(distance_ranges))

    best = []

    biggest_join = 1

    for n in range(len(data)): # iterate over all rank chunks
        best_result0, best_score = calc_score(join(data[:n+1]), distance_ranges)  # join all rank chunks into a single chunk
        best_result = {(0,n+1): best_result0}
        for last_cut in range(1, n+1): # last_cut is called K in the documentation
            curr_join = n + 1 - last_cut
            if curr_join > biggest_join * 2:
                continue
            partial_best_result, partial_best_score = best[last_cut-1]
            fwd_best_result, fwd_best_score = calc_score(join(data[last_cut:n+1]), distance_ranges)
            curr_score = partial_best_score + fwd_best_score
            #print("SC", n, last_cut, curr_score, partial_best_score, fwd_best_score, best_score)
            if curr_score > best_score:
                #print(n, last_cut, best_score)
                best_score = curr_score
                best_result = partial_best_result.copy()
                best_result[(last_cut, n+1)] = fwd_best_result
                if curr_join > biggest_join:
                    biggest_join = curr_join
        #print(dist_array_file, n, best_score)
        best.append((best_result, best_score))

    best_result, _ = best[-1]
    keys, values = [], []
    for k,v in best_result.items():
        keys.append([rank_ranges[k[0]], rank_ranges[k[1]]])
        values.append(v)
    potential = {
        "rank_chunks": keys,
        "distance_bins": values,
    }
    json.dump(potential, open(outfile, "w"), indent=2)

import multiprocessing
pool = multiprocessing.Pool(args.nparallel)
pool_args = []
for receptor_type in range(1, 99):
    for ligand_type in range(1, 99):
        f = "{}-{}.npy".format(receptor_type, ligand_type)
        dist_array_file = os.path.join(args.distance_array_dir, f)
        if os.path.exists(dist_array_file):
            f = "{}-{}.json".format(receptor_type, ligand_type)
            outfile = os.path.join(args.output_dir, f)
            pool_args.append((dist_array_file, outfile))
pool.starmap(make_histogram, pool_args)
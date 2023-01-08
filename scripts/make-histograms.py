"""Makes histograms from distance arrays.
Input: See build-distance-arrays.py for details
Second input: a JSON file containing the rank chunk list

One distance array must be present for every atomtype-atomtype combination
Each distance array contains small bins with (nearnative, nonnative) counts.
This is done separately for each chunk of ATTRACT ranks.

make-histograms first joins the rank chunks into bigger chunks. 
The rank chunk joinings are derived from a JSON input file 
In principle, all possible rank chunk joinings are tried
(with a speed-up, see below). 
For every rank chunk joining, a step potential consisting of distance bins (with a log-odds score for each bin)
is then automatically derived, by merging small bins until at least a minimum of near-natives are present in
the merged bin.

make-histograms-fixrankchunks should run within seconds.

"""
import os
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

p.add_argument("rank_chunks",
    help="A JSON file containing a list of rank chunks ( e.g [[0, 1000], [1000,2000]])"
)

p.add_argument("output_dir",
    help="The directory where to store the generated histograms, with names such as 1-32.json "
)

p.add_argument("--nparallel", type=int, default=None, help="Number of arrays to compute in parallel. By default, all cores are used.")

args = p.parse_args()

minimum_near_natives = args.minimum_near_natives

def get_bins(dist, minimum_near_natives):
    near_natives = dist[:, 0]
    non_natives = dist[:, 1]
    tot_near_natives = near_natives.sum()
    if not tot_near_natives:
        return []
    bins = []
    pos = 0
    while 1:
        nat, nonnat = 0, 0
        while 1:
            if pos == len(dist):
                break   
            nat += near_natives[pos]
            nonnat += non_natives[pos]
            pos += 1
            if nat >= minimum_near_natives:
                break
        bins.append(pos)
        if pos == len(dist):
            break
    return bins

def get_values(dist, distance_ranges, bins):
    assert len(dist) == len(distance_ranges)
    if not bins:
        return []
    distance_ranges = list(distance_ranges[1:]) + [999]
    near_natives = dist[:, 0]
    non_natives = dist[:, 1]
    tot_near_natives = near_natives.sum()
    tot_non_natives = non_natives.sum()
    if tot_near_natives == 0:
        baseline = -3
    elif tot_non_natives == 0:
        baseline = 3
    else:
        baseline = log(tot_near_natives) - log(non_natives.sum())    
    ranges, potential = [], []
    pos = 0
    for pos2 in bins:
        nat = near_natives[pos:pos2].sum()
        nonnat = non_natives[pos:pos2].sum() 
        if nat == 0:
            log_odds = -3
        elif nonnat == 0:
            log_odds = 3
        else:
            log_odds = log(nat) - log(nonnat) - baseline
        if log_odds < -3: log_odds = -3
        if log_odds > 3: log_odds = 3
        potential.append(log_odds)
        ranges.append(float(distance_ranges[pos2-1]))
        pos = pos2
    return list(zip(ranges, potential))

def make_histogram(dist_array_file, rank_chunks, outfile):
    dist_array = np.load(dist_array_file) # distance array, generated with build-distance-arrays.py
    distance_ranges = [0] + dist_array["distance_ranges"].tolist()
    rank_ranges = [0] + dist_array["rank_ranges"].astype(int).tolist()
    data = dist_array["data"].astype(np.uint32)  # counts
    assert data.shape[0] == len(rank_ranges)-1, (data.shape, len(rank_ranges)-1)
    assert data.shape[1] == len(distance_ranges), (data.shape, len(distance_ranges))

    bins = get_bins(data.sum(axis=0), minimum_near_natives)

    values = []
    for start, end in rank_chunks:
        try:
            startpos = rank_ranges.index(start)
        except ValueError:
            raise ValueError("{} not in rank ranges".format(start)) from None
        try:
            endpos = rank_ranges.index(end)
        except ValueError:
            if end == rank_chunks[-1][1]:
                endpos = len(rank_ranges) - 1
            else:
                raise ValueError("{} not in rank ranges".format(end)) from None
        currdata = data[startpos:endpos].sum(axis=0)
        currvalues = get_values(currdata, distance_ranges, bins)
        values.append(currvalues)

    potential = {
        "rank_chunks": rank_chunks,
        "distance_bins": values,
    }
    json.dump(potential, open(outfile, "w"), indent=2)

rank_chunks = json.load(open(args.rank_chunks))


import multiprocessing
if args.nparallel is None or args.nparallel > 1:
    pool = multiprocessing.Pool(args.nparallel)
pool_args = []
for receptor_type in range(1, 99):
    for ligand_type in range(1, 99):
        f = "{}-{}.npy".format(receptor_type, ligand_type)
        dist_array_file = os.path.join(args.distance_array_dir, f)
        if os.path.exists(dist_array_file):
            f = "{}-{}.json".format(receptor_type, ligand_type)
            outfile = os.path.join(args.output_dir, f)
            pool_args.append((dist_array_file, rank_chunks, outfile))

if args.nparallel is None or args.nparallel > 1:
    pool.starmap(make_histogram, pool_args)
else:
    for args in pool_args:
        make_histogram(*args)
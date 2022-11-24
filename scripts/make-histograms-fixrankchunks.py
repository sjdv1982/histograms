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
Therefore, there is a one-to-one relationship between step potential and joining
The best step potential is selected by discriminative power.
That power is computed as 

    abs(log_odds) * abs(nat-nat_expected), where:
        nat is the number of near-natives in the bin
        nat_expected is the expected number of near-natives in that bin.
          this is computed as fnat * (nat + nonnat)
          where fnat is sum(nat)/sum(nonnat) for all bins together
        log_odds is the log-odds score for that bin
    summed over all bins

In terms of speed up:
  The best joining of chunk 1..N is computed as the best joining of 1..K, plus a single chunk K..N,
  considering all values K between 1 and N. 
  The best joining of 1..K will have been computed before, as N is incremented from 1 to len(rank_chunks)
  This is just an optimization (no approximation = dynamic programming)

make-histograms usually runs within seconds.

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
p.add_argument("nparameters", type=int,
    help="The number of free parameters. This corresponds to bin boundaries. The final number of bins will be this + the number of rank chunks"
)

p.add_argument("rank_chunks",
    help="A JSON file containing a list of rank chunks ( e.g [[0, 1000], [1000,2000]])"
)

p.add_argument("output_dir",
    help="The directory where to store the generated histograms, with names such as 1-32.json "
)

p.add_argument("--nparallel", type=int, default=None, help="Number of arrays to compute in parallel. By default, all cores are used.")

args = p.parse_args()

nparameters = args.nparameters

def score(nat, nonnat, fnat):
    total = nat + nonnat
    if total == 0:
        return 0
    if nat > 0:
        log_odds = log((nat/total) / fnat)
    else:
        log_odds = -1000
    nat_expected = fnat * total
    return abs(log_odds) * abs(nat - nat_expected)

def optimize_boundaries(data, fnat):
    maxpar = min(nparameters, len(data)-1)
    result = []
    cumul = np.cumsum(data, axis=0)
    all_nat, all_nonnat = cumul[-1]
    result.append(((), score(*cumul[-1], fnat)))
    if maxpar == 0:
        return result
    best0 = [((),score(*c, fnat)) for c in cumul]    
    all_bestseries = [best0]
    last_bestseries = best0
    for i in range(1, maxpar+1):
        score_pre = last_bestseries[0][1]
        score_post = score(*data[i], fnat)
        curr_bestseries = [
            ((last_bestseries[0][0] + (i-1,)), score_pre + score_post)
        ]
        for n in range(i+2, len(data)+1):
            curr_best = None
            for nn in range(i+1, n+1):
                pos = nn - i
                score_pre = last_bestseries[pos][1]
                score_post = score(*data[nn:n].sum(axis=0), fnat)
                curr_best_score = score_pre + score_post 
                if curr_best is None or curr_best_score > curr_best[1]:
                    curr_best = (last_bestseries[pos][0] + (nn-1,), curr_best_score)
            curr_bestseries.append(curr_best)
        all_bestseries.append(curr_bestseries)
        last_bestseries = curr_bestseries
    return all_bestseries

def make_histogram(dist_array_file, rank_chunks, outfile):
    dist_array = np.load(dist_array_file) # distance array, generated with build-distance-arrays.py
    distance_ranges = [0] + dist_array["distance_ranges"].tolist()
    rank_ranges = [0] + dist_array["rank_ranges"].astype(int).tolist()
    data = dist_array["data"].astype(np.uint32)  # counts
    assert data.shape[0] == len(rank_ranges)-1, (data.shape, len(rank_ranges)-1)
    assert data.shape[1] == len(distance_ranges), (data.shape, len(distance_ranges))

    tot_near_natives = data[:, :, 0].sum()
    tot_non_natives = data[:, :, 1].sum()
    fnat = tot_near_natives / (tot_near_natives + tot_non_natives)

    data_chunks = []
    for start, end in rank_chunks:
        try:
            startpos = rank_ranges.index(start)
        except IndexError:
            raise ValueError("{} not in rank ranges".format(start)) from None
        try:
            endpos = rank_ranges.index(end)
        except IndexError:
            raise ValueError("{} not in rank ranges".format(end)) from None
        data_chunks.append(data[startpos:endpos].sum(axis=0))
    best_results = [optimize_boundaries(data_chunk, fnat) for data_chunk in data_chunks]

    
    best_chunk_scores = []
    for best_result in best_results:
        curr = []
        for k in range(nparameters+1):
            curr.append(best_result[k][-1][1])
        best_chunk_scores.append(curr)
    
    def select(tail_best_chunk_scores, remaining):
        if len(tail_best_chunk_scores) == 1:
            return tail_best_chunk_scores[0][remaining], [remaining]
        best = None
        for n in range(remaining):
            curr_score = tail_best_chunk_scores[0][n]
            tail = select(tail_best_chunk_scores[1:], remaining-n)
            tail_score = curr_score + tail[0]
            if best is None or tail_score > best[0]:
                best = tail_score, [n] + tail[1]
        return best

    best = select(best_chunk_scores,nparameters)
    best_nboundaries = best[1]
    best_boundaries = [best_results[b][k][-1][0] for b,k in enumerate(best_nboundaries)]
    
    print("FINAL", best_nboundaries, best_boundaries)
    bestscores = [best_results[b][k][-1][1] for b,k in enumerate(best_nboundaries)] 
    #print(best[0], bestscores, sum(bestscores))
    values = []
    for data_chunk, chunk_best_boundaries in zip(data_chunks, best_boundaries):
        start = 0
        #sc = 0            
        curr_values = []
        for boundary_pos in chunk_best_boundaries:
            nat, nonnat = data_chunk[start:boundary_pos+1].sum(axis=0)
            #sc0 = score(nat, nonnat, fnat)
            threshold = distance_ranges[boundary_pos + 1]
            if nat > 0:
                log_odds = log((nat/(nat + nonnat)) / fnat)
                log_odds += log(fnat)  # calibrate to get absolute log(p) of being native
            else:
                log_odds = -1000
            curr_values.append((threshold, log_odds))
            #sc += sc0
            start = boundary_pos+1
        #sc += score(*data_chunk[start:].sum(axis=0), fnat)
        #print(sc)
        if not chunk_best_boundaries:
            curr_values = [[0, log(fnat)]]
        values.append(curr_values)
    potential = {
        "rank_chunks": rank_chunks,
        "distance_bins": values,
    }
    json.dump(potential, open(outfile, "w"), indent=2)

rank_chunks = json.load(open(args.rank_chunks))


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
            pool_args.append((dist_array_file, rank_chunks, outfile))

pool.starmap(make_histogram, pool_args)
import sys
import numpy as np
from math import log

dist_array = np.load(sys.argv[1]) # distance array, generated with build-distance-arrays.py
minimum_near_natives = int(sys.argv[2]) # the minimum number of near-native contacts in each histogram bin
outfile = sys.argv[3] # JSON file

distance_ranges = [0] + dist_array["distance_ranges"].tolist()
rank_ranges = [0] + dist_array["rank_ranges"].tolist()
data = dist_array["data"]
assert data.shape[0] == len(rank_ranges)-1, (data.shape, len(rank_ranges)-1)
assert data.shape[1] == len(distance_ranges), (data.shape, len(distance_ranges))

def join(disc_dists):
    result = disc_dists[0].copy()
    result[0] = disc_dists[:, 0].sum()
    result[1] = disc_dists[:, 1].sum()
    return result

def calc_score(disc_dist):
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

best = []

biggest_join = 1

for n in range(len(data)):
    best_result0, best_score = calc_score(join(data[:n+1]))
    best_result = {(0,n+1): best_result0}
    for last_cut in range(1, n+1):
        curr_join = n + 1 - last_cut
        if curr_join > biggest_join * 2:
            continue
        partial_best_result, partial_best_score = best[last_cut-1]
        fwd_best_result, fwd_best_score = calc_score(join(data[last_cut:n+1]))
        curr_score = partial_best_score + fwd_best_score
        #print("SC", n, last_cut, curr_score, partial_best_score, fwd_best_score, best_score)
        if curr_score > best_score:
            #print(n, last_cut, best_score)
            best_score = curr_score
            best_result = partial_best_result.copy()
            best_result[(last_cut, n+1)] = fwd_best_result
            if curr_join > biggest_join:
                biggest_join = curr_join
    print(n, best_score)
    best.append((best_result, best_score))

import json
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
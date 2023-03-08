"""
220 cases
TO REMOVE: 1A9N-11-CCU (#7), 4QQB-11-CAC (#166), 4BS2-1-UGU (#131)
=> 217 cases

ATTRACT: 46 successes
  in the top 20pct (threshold: 60 % of all natives)
awk '$1 > 60' ATTRACT-5A-top20pct.txt  | wc -l

select 4: 96 successes
"""

import numpy as np
import os

THRESHOLD = 60   
PCT_THRESHOLD = 5 
error_runs = [7, 166,131]
pct = np.loadtxt("histo-5A-top5pct.txt").astype(float)

run_names = [l.strip() for l in open("curben.in-elem")]
scoresets = []
aliases = {}
for l in open("/home/sjoerd/histograms-data/rescore-div60-merged/aliases.txt"):
    ll = l.split()
    scoresets.append(ll[0])
    aliases[ll[0]] = ll[-1]

    assert len(run_names) == len(pct), (len(run_names), len(pct))

success = (pct >= THRESHOLD)
to_select = 4

def get_rank_natives(run1, run2, threshold):
    run1f = run1.replace(" ", "-")
    run2f = run2.replace(" ", "-")
    fname = os.path.join("/home/sjoerd/histograms-data/rescore-div60-merged/{}/{}.rank-natives".format(run1f, run2f))
    with open(fname) as f:
        lines = f.readlines()
    nstruc = int(lines[0].split()[-1])
    thresh = threshold/100 * nstruc
    lines2 = [l.split() for l in lines[1:]]
    indices = [int(l[0]) for l in lines2]
    ranks = [int(l[1]) for l in lines2]
    nnat = len(indices)
    return {ind for ind, r in zip(indices, ranks) if r <= thresh}, nnat



nsuccess = 0
done = set()
for redundancy_set in range(1000):
    f = "redundancy/redundancy-set-{}-indices.txt".format(redundancy_set+1)
    if not os.path.exists(f):
        break
    lines = open(f).readlines()
    test_set = [int(ll) for ll in lines[1].split()]
    training_set = [int(ll) for ll in lines[3].split()]
    rows_training_set = [n for n in range(len(pct)) if n+1 in training_set and n+1 not in error_runs]
    rows_test_set = [n for n in range(len(pct)) if n+1 in test_set and n+1 not in error_runs]
    selected =  np.arange(4, dtype=int) + 4 * redundancy_set
    success_submatrix = success[rows_test_set]
    success_testset_minimum = success_submatrix[:, selected].max(axis=1)    
    nsuccess_testset_minimum = success_testset_minimum.sum()

    nsuccess_trainset_minimum = success[rows_training_set][:, selected].max(axis=1).sum()
    nsuccess_trainset_minimum_CHECK = 0
    partitions = [[] for n in range(to_select)]
    for row in rows_training_set:
        p = pct[row, selected]
        if p.max() >= THRESHOLD:
            nsuccess_trainset_minimum_CHECK += 1
            partitions[p.argmax()].append(row)
    assert nsuccess_trainset_minimum_CHECK == nsuccess_trainset_minimum

    rows_test_set = rows_training_set ###
    nsuccess_testset = 0
    '''
    for column_nr, column in enumerate(selected):
        print(scoresets[column], aliases[scoresets[column]])
    '''
    for rownr, row in enumerate(rows_test_set):
        all_natives = set()
        nnat = None
        for column_nr, column in enumerate(selected):
            natives, nnat0 = get_rank_natives(run_names[row], aliases[scoresets[column]], PCT_THRESHOLD)
            if nnat is None:
                nnat = nnat0
            else:
                assert nnat == nnat0
            all_natives.update(natives)
        combined_pct = len(all_natives)/nnat * 100
        if combined_pct > THRESHOLD:
            nsuccess_testset += 1
    print(redundancy_set+1, nsuccess_testset, len(rows_test_set), "({:.2f})".format(100*nsuccess_testset/len(rows_test_set)))
    #print("*" * 50)
    #print()
    nsuccess += nsuccess_testset
    for run in test_set:
        done.add(run)

assert len(done) == len(pct)
print("Successes:", nsuccess)

'''
# Partition on the full set (for coneptual paper)
# 129/217 successes vs 96/217 for train/test

tokeep = [n for n in range(len(pct)) if n+1 not in error_runs]
success = success[tokeep][:, tokeep]
pct = pct[tokeep][:, tokeep]
run_names = [r for rnr, r in enumerate(run_names) if rnr in tokeep]

selected = sorted(get_selection(success))
success_minimum = success[:, selected].max(axis=1)    
nsuccess_minimum = success_minimum.sum()
pct2 = pct[:, selected]
pct2_sum = pct2.sum(axis=1)
nsuccess = 0
for row in range(len(success)):
    if success_minimum[row]:
        nsuccess += 1
        continue
    if pct2_sum[row] < THRESHOLD:
        continue
    all_natives = set()
    nnat = None
    for column_nr, column in enumerate(selected):
        if pct2[row, column_nr] == 0:
            continue
        natives, nnat0 = get_rank_natives(run_names[row], run_names[column], PCT_THRESHOLD)
        if nnat is None:
            nnat = nnat0
        else:
            assert nnat == nnat0
        all_natives.update(natives)
    combined_pct = len(all_natives)/nnat * 100
    #print("COMBINE", pct2[row], pct2_sum[row], combined_pct)
    if combined_pct > THRESHOLD:
        nsuccess += 1
print()
print("Success on full set (no train vs test):", nsuccess)
print()
print("# Selection")
for selnr, sel in enumerate(selected):
    print("#{}".format(selnr+1), run_names[sel])
print()

print("# Partition")
print("# complex fragnr frag best_partition pct_best_partition")
for n in range(len(success)):
    p = pct[n, selected]
    print(run_names[n], p.argmax()+1, p.max(), p)
'''
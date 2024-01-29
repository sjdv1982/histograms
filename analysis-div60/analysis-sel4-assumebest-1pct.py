"""
220 cases
TO REMOVE: 1A9N-11-CCU (#7), 4QQB-11-CAC (#166), 4BS2-1-UGU (#131)
=> 217 cases
"""

import numpy as np
import os

THRESHOLD = 60
PCT_THRESHOLD = 5
error_runs = [7, 166,131]
pct = np.loadtxt("histo-5A-top5pct.txt").astype(float)
pct_1pct = np.loadtxt("histo-5A-top1pct.txt").astype(float)
assert pct.shape[0] == pct.shape[1], pct.shape

run_names = [l.strip() for l in open("curben.in-elem")]
assert len(run_names) == len(pct), (len(run_names), len(pct))

success = (pct >= THRESHOLD)
to_select = 4

def get_rank_natives(run1, run2, threshold):
    run1f = run1.replace(" ", "-")
    run2f = run2.replace(" ", "-")
    fname = os.path.join("/home/sjoerd/histograms-data/rescore-div60/{}/{}.rank-natives".format(run1f, run2f))
    with open(fname) as f:
        lines = f.readlines()
    nstruc = int(lines[0].split()[-1])
    thresh = threshold/100 * nstruc
    lines2 = [l.split() for l in lines[1:]]
    indices = [int(l[0]) for l in lines2]
    ranks = [int(l[1]) for l in lines2]
    nnat = len(indices)
    return {ind for ind, r in zip(indices, ranks) if r <= thresh}, nnat


def select_columns(suc, score, best_score, selection, mincol, to_select):
    colsum = suc[:, mincol:].sum(axis=0)
    cols = np.argsort(-colsum)
    optimum = colsum[cols[:to_select]].sum()
    if score + optimum <= best_score:
        return None, None
    all_best_selection = None
    for col0 in cols:
        col = col0 + mincol
        next_mincol = col+1
        new_success = suc[:, col]
        next_score = score + new_success.sum()
        next_suc = suc.copy()
        next_suc[new_success, :] = 0
        next_suc[:, col] = 0
        next_selection = selection + [col]
        if to_select == 1:
            curr_best_score = next_score
            curr_selection = next_selection
        else:
            curr_selection, curr_best_score = select_columns(next_suc, next_score, best_score, next_selection, next_mincol, to_select-1)
        if curr_best_score is not None and curr_best_score > best_score:
            all_best_selection = curr_selection
            best_score = curr_best_score
    if all_best_selection is None:
        #if to_select > 6:
        #    print(to_select, best_score)
        return None, None
    else:
        #print(to_select, best_score, old_best_score)
        return all_best_selection, best_score

def get_selection(success_submatrix):
    colsum = success_submatrix.sum(axis=0)
    success2 = success_submatrix[:, np.argsort(-colsum)]
    selected0, best = select_columns(success2, 0, 0, [], 0, to_select)
    if best is None:
        selected = [0] * to_select
        best = 0
    else:
        selected = np.argsort(-colsum)[selected0]
    best2 = success_submatrix[:, selected].max(axis=1).sum()
    assert best == best2, (best, best2)
    return selected

stats = []
done = set()
for redundancy_set in range(1000):
    f = "redundancy/redundancy-set-{}-indices.txt".format(redundancy_set+1)
    if not os.path.exists(f):
        break
    lines = open(f).readlines()
    test_set = [int(ll) for ll in lines[1].split()]
    training_set = [int(ll) for ll in lines[3].split()]
    columns_to_keep = [n for n in range(len(pct)) if n+1 not in test_set]
    success_submatrix0 = success[:, columns_to_keep]
    rows_training_set = [n for n in range(len(pct)) if n+1 in training_set and n+1 not in error_runs]
    rows_test_set = [n for n in range(len(pct)) if n+1 in test_set and n+1 not in error_runs]
    selected0 = sorted(get_selection(success_submatrix0[rows_training_set]))
    ###selected0 = get_selection(success_submatrix0[rows_test_set]) #to test overfitting...
    if redundancy_set == 0:
        selected0 = [knr for knr,k in enumerate(columns_to_keep) \
            if k in [174,  63,  71, 150]]
    success_submatrix = success_submatrix0[rows_test_set]
    success_testset_minimum = success_submatrix[:, selected0].max(axis=1)    
    nsuccess_testset_minimum = success_testset_minimum.sum()
    pct_testset = pct[:, columns_to_keep][rows_test_set, :][:, selected0]
    pct_testset_1pct = pct_1pct[:, columns_to_keep][rows_test_set, :][:, selected0]
    pct_sum_testset = pct_testset.sum(axis=1)
    selected = np.array(columns_to_keep)[selected0]
    selected_run_names =  [rn for i,rn in enumerate(run_names) if i in selected]

    nsuccess_trainset_minimum = success_submatrix0[rows_training_set][:, selected0].max(axis=1).sum()
    nsuccess_trainset_minimum_CHECK = 0
    partitions = [[] for n in range(to_select)]
    for row in rows_training_set:
        p = pct[row, selected]
        if p.max() >= THRESHOLD:
            nsuccess_trainset_minimum_CHECK += 1
            partitions[p.argmax()].append(row)
    assert nsuccess_trainset_minimum_CHECK == nsuccess_trainset_minimum

    for rownr, row in enumerate(rows_test_set):
        assert len(pct_testset[rownr]) == to_select
        x = pct_testset_1pct[rownr].argmax()
        print(run_names[row], pct_testset_1pct[rownr].max(), selected_run_names[x])
        stats.append((run_names[row], pct_testset_1pct[rownr].max()))

    for run in test_set:
        done.add(run)

assert len(done) == len(pct)

cmpl = {}
for run_name, p in stats:
    c = run_name.split()[0]
    if p > cmpl.get(c, -1):
        cmpl[c] = p

for c in cmpl:
    print(c, cmpl[c])
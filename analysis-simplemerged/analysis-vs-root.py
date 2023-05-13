# comparison with pure non-merged

THRESHOLD = 60   
PCT_THRESHOLD = 5 

import os
import numpy as np

roots = [45, 151, 175, 196]


def get_rank_natives_OLD(run1, run2, threshold):
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

def get_rank_natives(run1, run2, threshold):
    run1f = run1.replace(" ", "-")
    run2f = run2.replace(" ", "-")
    fname = os.path.join("/home/sjoerd/histograms-data/rescore-simplemerge/{}/{}.rank-natives".format(run1f, run2f))
    with open(fname) as f:
        lines = f.readlines()
    nstruc = int(lines[0].split()[-1])
    thresh = threshold/100 * nstruc
    lines2 = [l.split() for l in lines[1:]]
    indices = [int(l[0]) for l in lines2]
    ranks = [int(l[1]) for l in lines2]
    nnat = len(indices)
    return {ind for ind, r in zip(indices, ranks) if r <= thresh}, nnat

run_names = [l.strip() for l in open("curben.in-elem")]

pat_new = "weight25-partition"

ncases = 0
nsuccess_old = 0
nsuccess_new = 0

for run_nr, run_name in enumerate(run_names):
    results_old, results_new = [], []
    for part in 1,2,3,4:
        try:
            pat_old = run_names[roots[part-1]-1]
            result_old = get_rank_natives_OLD(run_name, pat_old, PCT_THRESHOLD)
            result_new = get_rank_natives(run_name, pat_new + str(part), PCT_THRESHOLD)
            assert result_old[1] == result_new[1]
            results_old.append(result_old)
            results_new.append(result_new)
            assert result_old[1] == results_old[0][1]
        except Exception:
            print("SKIP", run_name)
            break
    else:
        nnat = results_new[0][1]
        all_natives_old = set()
        for r in results_old:
            all_natives_old.update(r[0])
        combined_pct_old = len(all_natives_old)/nnat * 100
        all_natives = set()
        for r in results_new:
            all_natives.update(r[0])
        combined_pct = len(all_natives)/nnat * 100

        print("{} {:.1f} {:.1f}".format(run_name, combined_pct, combined_pct_old))
        ncases += 1
        if combined_pct_old >= THRESHOLD:
            nsuccess_old += 1        
        if combined_pct >= THRESHOLD:
            nsuccess_new += 1
        

print(ncases, nsuccess_new, nsuccess_old)            
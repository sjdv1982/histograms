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
assert pct.shape[0] == pct.shape[1], pct.shape

run_names = [l.strip() for l in open("curben.in-elem")]
assert len(run_names) == len(pct), (len(run_names), len(pct))
cases = [rnr for rnr, r in enumerate(run_names) if rnr+1 not in error_runs]
case_names = [r for rnr, r in enumerate(run_names) if rnr+1 not in error_runs]
pct = pct[cases]

selected_names = ['1M5K 4 GCA', '4N0T 10 AGA', '5MPG 1 UAG', '6DCL 5 UUA']
selected = [run_names.index(s) for s in selected_names]

stats = []
for case, case_name in enumerate(case_names):
    mx = -1
    for sel in selected:
        cur = pct[case, sel]
        if cur > mx:
            mx = cur
    print(case_name, mx)
    stats.append((case_name, mx))

cmpl = {}
for run_name, p in stats:
    c = run_name.split()[0]
    if p > cmpl.get(c, -1):
        cmpl[c] = p

for c in cmpl:
    print(c, cmpl[c])
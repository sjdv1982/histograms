to_remove = [7, 166,131]
import os
runs = [l.strip() for l in open("curben.in-elem")]
f1 = open("ATTRACT-5A-top5pct.txt", "w")

for runnr, run in enumerate(runs):
    if runnr+1 in to_remove:
        continue
    runf = run.replace(" ", "-")
    for run2 in [run]:
        run2f = run2.replace(" ", "-")
        fname = os.path.join("/home/sjoerd/histograms-data/rescore-div60/{}/{}.rank-natives".format(runf, run2f))
        if runnr+1 in to_remove:
            pct = -1
            pct_attract = -1
        else:
            assert os.path.exists(fname), fname
            with open(fname) as f:
                lines = f.readlines()
            nstruc = int(lines[0].split()[-1])
            thresh_histo = 0.05 * nstruc
            thresh_attract = 0.05 * nstruc 
            lines2 = [l.split() for l in lines[1:]]
            nnat = len(lines2)
            ranks_attract = [int(l[0]) for l in lines2]
            pct_attract = len([r for r in ranks_attract if r <= thresh_attract]) / nnat * 100
            pct_attract = "%.1f" % pct_attract
    print(run, pct_attract, file=f1)
    f1.flush()

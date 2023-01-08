import numpy as np
import sys
scorefile=sys.argv[1]
complex, frag = scorefile.split("-")[:2]
rmsds=np.loadtxt("../{}/rmsd/frag{}.lrmsd".format(complex, frag), np.float32)[:, 1]
scores=np.loadtxt(scorefile,np.float32)
ranks = np.argsort(-scores)
print("# histograms vs ATTRACT")
for rmsd in 2, 3, 5:
    print("# RMSD: {}".format(rmsd))
    nat=(rmsds<rmsd)
    print("natives: {:d}".format(nat.sum()))
    for top in 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000:
        top_histo = nat[ranks[:top]].sum()
        top_attract = nat[:top].sum()
        print(rmsd, top, top_histo, top_attract)
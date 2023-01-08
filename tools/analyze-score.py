import numpy as np
import sys
rank_native_file=sys.argv[1]
rank_natives = np.loadtxt(rank_native_file)
ranks_ATTRACT = rank_natives[:, 0].astype(int)
ranks_histo = rank_natives[:, 1].astype(int)
rmsds = rank_natives[:, 2]
print("# histograms vs ATTRACT")
for rmsd in 2, 3, 5:
    print("# RMSD: {}".format(rmsd))
    nat=(rmsds<rmsd)
    ranks_ATTRACT2 = ranks_ATTRACT[nat]
    ranks_histo2 = ranks_histo[nat]
    print("natives: {:d}".format(nat.sum()))
    for top in 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000:
        top_histo = (ranks_histo2 <= top).sum()
        top_attract = (ranks_ATTRACT2 <= top).sum()
        print(rmsd, top, top_histo, top_attract)

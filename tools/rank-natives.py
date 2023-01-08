import numpy as np
import sys
scorefile=sys.argv[1]
complex, frag = scorefile.split("-")[:2]
rmsds=np.loadtxt("../{}/rmsd/frag{}.lrmsd".format(complex, frag), np.float32)[:, 1]
scores=np.loadtxt(scorefile,np.float32)
assert len(rmsds) == len(scores)
ranks= (-scores).argsort().argsort()
print("# RMSD ranks < 5, ATTRACT vs reranked, nstruc: {}".format(len(scores)))
rmsds_ranksorted = rmsds[ranks]
natives = np.where(rmsds < 5)[0]
for native in natives:
    print(native + 1, ranks[native]+1, "%.3f" % rmsds[native])

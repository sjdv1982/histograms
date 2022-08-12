import numpy as np, argparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("file",
help="""Name of the distance array file""")

p.add_argument("--rank-chunk-index", dest="rchunk", type=int, help= """Rank chunk to print out. If not supplied, print general information instead""")

args = p.parse_args()

dist_array = np.load(args.file) # distance array, generated with build-distance-arrays.py
distance_ranges = dist_array["distance_ranges"]
rank_ranges = dist_array["rank_ranges"]
if args.rchunk is None:
    print("Distance ranges:", distance_ranges)
    print("Rank ranges:", rank_ranges)
    exit(0)
chunk = dist_array["data"][args.rchunk-1]
distance_ranges = [0] + distance_ranges.tolist()
print("distance- distance+ #natives #nonnatives")
for n in range(len(chunk)):
    vmin = distance_ranges[n]
    vmax = 99
    if n < len(chunk)-1:
        vmax = distance_ranges[n+1]
    d = chunk[n]
    print("%.3f %.3f %d %d" % (vmin, vmax, d[0], d[1]))
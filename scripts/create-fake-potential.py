# reads in an idealized rank-chunks.json, an existing histogram (e.g. 1-32.json),
# and writes out a histogram with the potential (distance_bins) of the first rank chunk copied into all others
import sys, json
histo = json.load(open(sys.argv[1]))
rank_chunks = json.load(open(sys.argv[2]))
first_distance_bins = histo["distance_bins"][0]
result = {
    "rank_chunks": rank_chunks,
    "distance_bins": [first_distance_bins] * len(rank_chunks),
}
json.dump(result, sys.stdout, indent=2)

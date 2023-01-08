
# rebases an existing set of histograms to use a fixed rank chunk scheme

import os
import numpy as np
import argparse
import json
import glob

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

p.add_argument("input_dir",
    help="The directory with the original histograms, with names such as 1-32.json "
)

p.add_argument("rank_chunks",
    help="A JSON file containing a list of rank chunks ( e.g [[0, 1000], [1000,2000]])"
)

p.add_argument("output_dir",
    help="The directory where to store the generated histograms, with names such as 1-32.json "
)

args = p.parse_args()

rank_chunks = json.load(open(args.rank_chunks))

histofiles = glob.glob(os.path.join(args.input_dir, "[0-9]*-[0-9]*.json"))
for histofile in histofiles:
    with open(histofile) as f:
        histo = json.load(f)
    distance_bins = []
    for r1 in rank_chunks:
        max_overlap = 0
        for r2nr, r2 in enumerate(histo["rank_chunks"]):
            borders = r1[0], r1[1], r2[0], r2[1]
            argsort = np.argsort(borders).tolist()
            if argsort not in ([0,1,2,3], [2,3,0,1]):
                # overlap between source and target
                overlap = borders[argsort[2]] - borders[argsort[1]]
                if overlap > max_overlap:
                    max_overlap = overlap
                    max_r2nr = r2nr
        if max_overlap > 0:
            distance_bins.append(histo["distance_bins"][max_r2nr])
        else:
            raise Exception(histofile, r1)
    potential = {
        "rank_chunks": rank_chunks,
        "distance_bins": distance_bins,
    }
    outfile = os.path.join(args.output_dir, os.path.split(histofile)[1])
    json.dump(potential, open(outfile, "w"), indent=2)
    print(outfile)
    

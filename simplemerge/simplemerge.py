# Training set #20, partition #1
root1 = 45
members1 = "# Members: 4 5 6 18 24 25 26 27 31 32 33 44 45 46 47 48 49 52 53 54 84 100 115 118 137 171 177 200 201 202 203 204"
members1 = [int(l) for l in members1.split()[2:]]


# Training set #20, partition #2
root2 = 151
members2 = "# Members: 41 58 65 68 69 70 78 79 85 86 110 111 113 132 133 134 136 151 152 153 183 184 185 190 211"
members2 = [int(l) for l in members2.split()[2:]]

# Training set #20, partition #3
root3 = 175
members3 = "# Members: 13 17 40 50 59 60 63 64 73 74 75 82 83 93 94 104 105 106 107 114 126 127 138 139 144 161 162 173 174 175 176 179 182 186 188 192 193 194 205"
members3 = [int(l) for l in members3.split()[2:]]

# Training set #20, partition #4
root4 = 196
members4 = "# Members: 1 2 3 62 66 67 72 89 90 91 119 121 122 135 156 172 180 195 196 197 198"
members4 = [int(l) for l in members4.split()[2:]]

roots = [root1, root2, root3, root4]
allmembers = [members1, members2, members3, members4]

import json, os
import numpy as np

cases = [l.split() for l in open("curben.in-elem")]
cases = [c for c in cases if len(c) == 3]

FIRST = None
histos = {}
d = "/home/sjoerd/histograms-data/histograms-rebased"
for casenr, (c, frag, seq) in enumerate(cases):
    dd = os.path.join(d, c, "{}-{}".format(frag, seq))
    if not os.path.exists(dd):
        print("WARNING: {} {} does not exist".format(casenr+1, dd))
        continue
    h = {}
    histos[casenr+1] = h
    for n in range(1, 31+1):
        for nn in range(32,48+1):
            f = "{}/{}-{}.json".format(dd, n, nn)
            if not os.path.exists(f):
                continue  
            with open(f) as ff:
                data = json.load(ff)
            if FIRST is None:
                FIRST = data
            assert data["rank_chunks"] == FIRST["rank_chunks"]
            assert len(data["distance_bins"]) == len(FIRST["rank_chunks"])
            h[n,nn] = [np.array(dat) for dat in data["distance_bins"]]

nchunks = len(FIRST["rank_chunks"])

allh = []
for h in histos.values():
    for hh in h.values():
        for hhh in hh:
            if len(hhh):
                allh.append(hhh)
distances = np.unique(np.concatenate(allh)[:, 0])

histos_disc = {}
for key, h in histos.items():
    hdisc = {}
    for key2, hh in h.items():
        hhdisc = []
        for hhh in hh:
            if not len(hhh):
                hhdisc.append(None)
                continue
            hhhdisc = np.empty(len(distances))
            hhhpos = np.searchsorted(hhh[:, 0], distances)
            try:
                hhhdisc = hhh[:, 1][hhhpos]
            except IndexError:
                vals = np.append(hhh[:, 1], 0)
                hhhdisc = vals[hhhpos]
            hhdisc.append(hhhdisc)
        hdisc[key2] = hhdisc
    histos_disc[key] = hdisc

for weight in 25, 50, 75, 100:
    w = weight / 100
    wmin = 1 - w
    for part in 1,2,3,4:
        d = "/home/sjoerd/histograms-data/histograms-simple/weight{}-partition{}".format(weight, part)
        os.makedirs(d, exist_ok=True)
        root = roots[part-1]
        members = allmembers[part-1]
        for n in range(1, 31+1):
            for nn in range(32,48+1):
                distance_bins = []
                skip = False
                for chunk in range(nchunks):
                    key = n, nn
                    merged = np.zeros(len(distances))
                    count = 0
                    for member in members:
                        if member not in histos_disc:
                            continue
                        h = histos_disc[member]
                        if key not in h:
                            continue
                        if h[key][chunk] is None:
                            continue
                        merged += h[key][chunk]
                        count += 1                    
                    h = histos_disc[root]
                    if key in h and wmin > 0:
                        merged *= (w / count)
                        merged += (wmin * h[key][chunk])
                    elif count == 0:                
                        skip = True
                        break
                    else:
                        merged /= count
                    mergedlist = []
                    for pos, dis in enumerate(distances):
                        if pos + 1 == len(distances) or merged[pos+1] != merged[pos]:
                            mergedlist.append((dis, merged[pos]))                         
                    distance_bins.append(mergedlist)
                if skip:
                    print("SKIP", weight, part, n, nn)
                    continue
                dic = {
                    "rank_chunks": FIRST["rank_chunks"],
                    "distance_bins": distance_bins,
                }
                f = "{}/{}-{}.json".format(d, n, nn)
                with open(f, "w") as ff:                    
                    json.dump(dic, ff)
                #print(f)
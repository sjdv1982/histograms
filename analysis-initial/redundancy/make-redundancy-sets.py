import os

seqs = []
code_to_aln = {}
for lnr, l in enumerate(open("my_rrm_alignment.fasta")):
    if l.startswith(">"): 
        code = l.split("_")[2]
        if code not in code_to_aln:
            code_to_aln[code] = []
        code_to_aln[code].append(len(seqs) + 1) # count start at 1
        continue
    l = l.strip()
    if len(l) == 2: continue # error in file
    seqs.append(l.upper())
    if len(l) != 243: print(l, len(l), lnr)


def get_seqid(seq1, seq2):
    assert len(seq1) == len(seq2)
    count = 0
    tot = 0
    for s1, s2 in zip(seq1, seq2):
        if s1 == "-" and s2 == "-":
            continue
        tot += 1
        if s1 == s2:
            count += 1
    return count/tot * 100

redundant = {}
for code in code_to_aln:
    redundant[code] = set()
    for code2 in code_to_aln:
        highseqid = False
        for s1 in code_to_aln[code]:
            for s2 in code_to_aln[code2]:
                seqid = get_seqid(seqs[s1-1], seqs[s2-1])
                if seqid > 40:
                    #print(code, code2, seqid)
                    highseqid = True
                    break
            if highseqid:
                break
        if highseqid:
            redundant[code].add(code2)


bestgroups = None

import random
random.seed(0)

#for n in range(10000): # doesn't change...
for n in range(1): # doesn't change...    
    done = set()
    groups = []
    order = list(code_to_aln.keys())
    random.shuffle(order)
    for code in order:
        if code in done:
            continue
        redund = redundant[code]
        testset = set()
        testset.add(code)
        trainset = {c for c in code_to_aln if c not in redund}
        for code2 in redund:
            if code2 in done:
                continue
            for code3 in redundant[code2]:
                if code3 in trainset:
                    break
            else:
                testset.add(code2)
        done.update(testset)
        groups.append((testset, trainset))

    if bestgroups is not None and len(bestgroups) <= len(groups):
        continue

    for group in groups:
        testset, trainset = group
        for code in testset:
            for code2 in redundant[code]:
                assert code2 not in trainset
    bestgroups = groups

groups = bestgroups
print("Redundancy sets: {}".format(len(groups)))

code_to_docking = {code: [] for code in code_to_aln}
for nr, l in enumerate(open("curben.in-elem")):
    code = l.split()[0]
    try:
        code_to_docking[code]
    except KeyError:
        print("{} not in alignment, consider it as non-redundant (manually verified for 3MOJ,3RW6)".format(code))
        testset = {code}
        trainset = {c for c in code_to_aln if c != code}
        for group in groups:
            group[1].add(code)
        groups.append((testset, trainset))
        code_to_docking[code] = []
    code_to_docking[code].append(nr+1)

os.system("rm -f redundancy-set-*.txt")
print("redundancy-set-*-codes.txt refer to PDB codes")
print("redundancy-set-*-indices.txt refer to indices (starting at 1) in curben.in-elem")
groups = sorted(groups, key=lambda group: -len(group[0]))
for groupnr, group in enumerate(groups):
    with open("redundancy-set-{}-codes.txt".format(groupnr+1), "w") as f:
        print("Test set", file=f)
        for code in sorted(group[0]):
            print(code,end=" ",file=f)
        print(file=f)
        print("Training set", file=f)
        for code in sorted(group[1]):
            print(code,end=" ",file=f)
        print(file=f)        
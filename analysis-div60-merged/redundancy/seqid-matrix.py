seqs = []
for lnr, l in enumerate(open("my_rrm_alignment.fasta")):
    if l.startswith(">"): continue
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

for i in range(len(seqs)):
    for j in range(len(seqs)):
        seqid = get_seqid(seqs[i], seqs[j])
        print(i+1, j+1, "{:.1f}".format(seqid))

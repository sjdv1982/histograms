import sys, os
import numpy as np
pdbfile = sys.argv[1]  # ligand PDB with atom types, e.g. library/GCA/conf-1.pdb
datfile = sys.argv[2]  # e.g. GCA.dat
liblist = sys.argv[3]  # e.g. library/GCA.list
output_pattern = sys.argv[4]  # e.g. coordinates/GCA
                              # Will create coordinates/GCA-32.npy, coordinates/GCA-33.npy, ...
                              # where 32,33 are the atom types

atomtypes = np.array([int(l[57:59]) for l in open(pdbfile) if l.startswith("ATOM")])
for atomtype in np.unique(atomtypes):
    indices = np.where(atomtypes == atomtype)[0] + 1
    inds = " ".join([str(i) for i in indices])
    outfile = output_pattern + "-" + str(atomtype) + ".npy"
    cmd = "python2 $ATTRACTDIR/dump_coordinates.py %s %s %s %d %s --ens 2 %s" % (datfile, pdbfile, outfile, len(indices), inds, liblist)
    print(cmd)
    os.system(cmd)
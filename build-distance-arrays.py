"""Builds Numpy arrays of near-native/non-native distance statistics.

Input is a complex directory.
It is assumed that the following files exist:
- <complex>/receptor-reduced.pdb . Receptor PDB. All docking poses are assumed to have the same receptor conformation
- <complex>/<motif>.rmsd . For each docked motif, the (ligand) RMSD file
- <complex>/coordinates/<motif>-<atomtype>.npy. For each docked motif, all the coordinates for that atom type.

All docking poses are assumed to be ranked by ATTRACT energy.
The number of poses must always be the same for each motif.

The output is one Numpy structure for each atomtype-atomtype pair
The structure contains a distance range array, a rank range array, and a data array

The distances are discretized at specified distance discretizations,
until specified maximum distances. Starting at zero, this leads to
N distance ranges, where the last distance range has infinity as the upper limit.
The N-1 upper limits are stored in the distance range array.

The docking poses are divided into ATTRACT rank chunks of a specified length.
K is the total number of chunks. The K upper rank limits of each chunk are
stored in the rank range array.

The shape of each data array is K * N bins of 2 integers.
Each bin contains the number of near-natives and the number of non-natives.
For each rank chunk, there is one bin for each distance discretization.

Numpy arrays are saved in distance-arrays/<receptor-atomtype>-<ligand-atomtype>.npy
"""
import sys, os, argparse
import numpy as np
import itertools

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("complex",
help="""Name of the complex directory.
It is assumed that the following files exist:
- <complex>/receptor-reduced.pdb . Receptor PDB. All docking poses are assumed to have the same receptor conformation
- <complex>/<motif>.rmsd . For each docked motif, the (ligand) RMSD file
- <complex>/coordinates/<motif>-<atomtype>.npy. For each docked motif, all the coordinates for that atom type.
""")

p.add_argument("--distances",type=float, required=True, nargs="*", help= """Distance thresholds for each discretization regime.
Must be equal to the number of discretization values.
For example, "--distances 1 7.5 25" specified three discretization regimes: from 0 to 1, from 1 to 7.5, and from 7.5 to 25 A.
Each of these regimes will have its own discretization value, leading to a certain number of bins.
The final bin is for distances beyond 25 A.
""")

p.add_argument("--discretizations",type=float, required=True, nargs="*", help= """Discretization values for each discretization regime.
Must be equal to the number of distance thresholds.
For example, "--discretizations 0.1 0.02 0.1 0.2" specifies values for four discretization regimes.
A discretization regime with range of 0 to 1 A and discretization value 0.1 A will have 10 bins.
Note that the final bin is always for distances beyond the last distance threshold.
""")

p.add_argument("--rank-chunk", type=int, required=True, help="""ATTRACT rank chunk size.
The docking poses are divided into chunks of this size.
Note that it is assumed that the docking poses are already sorted by ATTRACT rank.""")

p.add_argument("--rmsd-thresholds", type=float, required=True, nargs=2, help="""RMSD tresholds.
The first threshold is the maximum RMSD for near-native structures
The second threshold is the minimum RMSD for non-native structures
Structures with an RMSD in between are ignored.
""")

p.add_argument("--nparallel", type=int, default=None, help="Number of arrays to compute in parallel. By default, all cores are used.")

args = p.parse_args()
assert len(args.distances) == len(args.discretizations)

bins = [0]
lower_threshold = 0
for upper_threshold, discretization in zip(args.distances, args.discretizations):
    threshold = lower_threshold
    while 1:
        threshold += discretization
        if threshold > upper_threshold:
            threshold = upper_threshold
        bins.append(threshold)
        if threshold == upper_threshold:
            break
    lower_threshold = upper_threshold
bins = np.array(bins, dtype=np.float32)

receptor = {}
for l in open(args.complex + "/receptor-reduced.pdb"):
    if not l.startswith("ATOM"):
        continue
    atomtype = int(l[57:59])
    x, y, z = float(l[30:38]),float(l[38:46]),float(l[46:54])
    if atomtype not in receptor:
        receptor[atomtype] = []
    receptor[atomtype].append((x,y,z))
for atomtype in receptor:
    receptor[atomtype] = np.array(receptor[atomtype])

all_receptor_atomtypes = set()
for aa in receptor:
    if aa == 99:
        continue
    all_receptor_atomtypes.add(aa)

bases = "A", "C", "G", "U"
all_motifs = ["".join(comb) for comb in itertools.product(bases, bases, bases)]

def rmsd_file(complex, motif):
    return "%s/%s.rmsd" % (complex, motif)

def coordinate_file(complex, motif, atomtype):
    return "%s/coordinates/%s-%s.npy" % (complex, motif, atomtype)

motifs = []
for motif in all_motifs:
    if os.path.exists(rmsd_file(args.complex, motif)):
        motifs.append(motif)

if not len(motifs):
    raise Exception("No docking file found for any motif")

nstruc = None
for receptor_atomtype in all_receptor_atomtypes:
    for ligand_atomtype in range(1, 99):
        for motif in motifs:
            cfile = coordinate_file(args.complex, motif, ligand_atomtype)
            if not os.path.exists(cfile):
                continue
            ligand_coordinates = np.load(cfile)
            nstruc = len(ligand_coordinates)
            break
        if nstruc is not None:
            break
    if nstruc is not None:
        break

if nstruc is None:
    raise Exception("No coordinate file found")

rank_chunks = list(range(0, nstruc, args.rank_chunk))
if rank_chunks[-1] != nstruc:
    rank_chunks.append(nstruc)

struc_type = np.dtype([
    ("distance_ranges", np.float32, len(bins)-1),
    ("rank_ranges", np.uint32, len(rank_chunks)-1),
    ("data", np.uint32, (len(rank_chunks)-1, len(bins), 2)),
])

def calculate(receptor_atomtype, ligand_atomtype):
    result = np.zeros(1, struc_type)[0]
    result["distance_ranges"] = bins[1:]
    result["rank_ranges"] = rank_chunks[1:]
    receptor_coordinates = receptor[receptor_atomtype]
    exists = False
    for motif in motifs:
        cfile = coordinate_file(args.complex, motif, ligand_atomtype)
        if not os.path.exists(cfile):
            continue
        exists = True
        rmsds = np.array([float(l) for l in open(rmsd_file(args.complex, motif)) if len(l.strip())])
        native_mask = (rmsds <= args.rmsd_thresholds[0])
        nonnative_mask = (rmsds >= args.rmsd_thresholds[1])
        del rmsds
        ligand_coordinates = np.load(cfile)
        assert len(ligand_coordinates) == nstruc, (cfile, len(ligand_coordinates), nstruc)
        for k in range(len(rank_chunks)-1):
            min_rank, max_rank = rank_chunks[k:k+2]
            curr_coordinates = ligand_coordinates[min_rank:max_rank]
            d = result["data"][k]
            for nat, nativeness_mask in enumerate((native_mask, nonnative_mask)):
                coordinates = curr_coordinates[nativeness_mask[min_rank:max_rank]]
                for coor in receptor_coordinates:
                    if len(coordinates):
                        dist = np.linalg.norm(coordinates - coor, axis=1)
                        curr_counts, _ = np.histogram(dist, bins, density=False)
                        d[:-1, nat] += curr_counts.astype(np.uint32)
                        d[-1, nat] += (dist > bins[-1]).sum().astype(np.uint32)
        print(receptor_atomtype, ligand_atomtype, motif)

    if exists:
        np.save("%s/distance-arrays/%s-%s.npy" % (args.complex, receptor_atomtype, ligand_atomtype), result)

import multiprocessing
pool = multiprocessing.Pool(args.nparallel)
args = []
for receptor_atomtype in all_receptor_atomtypes:
    for ligand_atomtype in range(1, 99):
        #calculate(receptor_atomtype, ligand_atomtype)
        args.append((receptor_atomtype, ligand_atomtype))
pool.starmap(calculate, args)
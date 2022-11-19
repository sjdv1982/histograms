# Currently not used, but can be used to write discretized coordinates to disk.
# The "discretize_coordinates" function has been copied into score-with-histograms-discrete.py and optimized
import sys
import os
import json
import numpy as np
from hashlib import sha3_256

def discretize_coordinates(coordinates, rank_chunks):
    bits10 = np.uint32(2**10-1)

    coor_min = np.floor(coordinates.min(axis=0).min(axis=0))
    coordinates -= coor_min
    rank_start = 0

    result_dtype = np.dtype([("index", np.uint32),("weight", np.uint8)],align=False)
    digit_indices = np.empty(coordinates.shape[:2] + (8,), dtype=result_dtype)

    all_digit_coors = []
    for rank_start, rank_end in rank_chunks:
        chunk = coordinates[rank_start:rank_end]    
        if not len(chunk):
            continue
        discrete = []
        for xyz in range(3):
            coor = chunk[:, :, xyz]
            minus = np.floor(coor).astype(np.uint32)
            plus = minus + 1
            assert plus.max() < 1024
            wminus = plus - coor
            wplus = coor - minus
            discrete.append(((minus, wminus),(plus, wplus)))
        # digit_chunk = np.empty(chunk.shape, np.uint32)
        digit_hash = np.empty((8,) + chunk.shape[:2], np.uint32)
        weights = np.empty((8,) + chunk.shape[:2])
        ind_xyz = 0
        for dcoorx, weightx in discrete[0]:
            # digit_chunk[:, :, 0] = dcoorx
            hash_x = (dcoorx << 20)
            for dcoory, weighty in discrete[1]:
                # digit_chunk[:, :, 1] = dcoory
                hash_y = (dcoory << 10)
                weightxy = weightx * weighty
                for dcoorz, weightz in discrete[2]:
                    # digit_chunk[:, :, 2] = dcoorz
                    digit_hash[ind_xyz] = hash_x + hash_y + dcoorz
                    weights[ind_xyz] = weightxy * weightz
                    ind_xyz += 1
        uniq, digit_index = np.unique(digit_hash, return_inverse=True)
        digit_coors = np.empty((len(uniq), 3))
        digit_coors[:, 0] = (uniq >> 20) & bits10
        digit_coors[:, 1] = (uniq >> 10) & bits10
        digit_coors[:, 2] = uniq & bits10
        digit_index = digit_index.reshape(digit_hash.shape).astype(np.uint32)
        print("Number of atoms: {}, unique on grid: {}".format(chunk.shape[0] * chunk.shape[1], len(uniq)), file=sys.stderr)
        digit_indices["index"][rank_start:rank_end] = np.moveaxis(digit_index, 0, -1)
        weights = np.round(weights * 255).astype(np.uint8)
        digit_indices["weight"][rank_start:rank_end] = np.moveaxis(weights, 0, -1)
        all_digit_coors.append(digit_coors)
    
    digit_coors = np.concatenate(all_digit_coors) + coor_min
    return digit_coors, digit_indices

def main():
    coordinates = sys.argv[1] # e.g GUG-32.npy
    ligand_atomtype = int(sys.argv[2]) # e.g. 32
    histogram_file = sys.argv[3] # e.g. histograms/6-38.json
    output_dir = sys.argv[4] # will write discrete-coordinates-X-Y.json, where:
    #  X contains the ligand type, and refers to the coordinates file
    #  Y is the checksum of the rank chunks inside the histogram file.
    #  If the output file already exists, nothing needs to be done. 

    coordinates = np.load(coordinates)
    nstruc = len(coordinates)

    with open(histogram_file) as f:
        histogram = json.load(f)

    rank_chunks0 = histogram["rank_chunks"]
    rank_chunks = []
    for rank_start, rank_end in rank_chunks0:
        rank_start, rank_end = int(rank_start), int(rank_end)
        if rank_start >= nstruc:
            continue
        if rank_end > nstruc:
            rank_end = nstruc
        rank_chunks.append((rank_start, rank_end))
    rankchunk_checksum = sha3_256(json.dumps(rank_chunks, indent=2).encode()).digest().hex()

    output_pattern = os.path.join(output_dir, "discrete-{}-{}".format(ligand_atomtype, rankchunk_checksum))
    outputfile_1 = output_pattern + "-coordinates.npy"
    outputfile_2 = output_pattern + "-indices.npy"
    for outputfile in (outputfile_1, outputfile_2):
        if os.path.exists(outputfile):
            print("{} already exists, nothing to do".format(outputfile))
            sys.exit(0)

    output_pattern = os.path.join(output_dir, "discrete-{}-{}".format(ligand_atomtype, rankchunk_checksum))
    outputfile_1 = output_pattern + "-coordinates.npy"
    outputfile_2 = output_pattern + "-indices.npy"
    

    digit_coors, digit_indices = discretize_coordinates(coordinates, rank_chunks)
    np.save(outputfile_1, digit_coors)
    np.save(outputfile_2, digit_indices)

if __name__ == "__main__":
    main()
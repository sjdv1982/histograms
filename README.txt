See in each script details about what inputs are needed.

1. calculate all per-fragment RMSD files in parallel (lrmsd-all.py)
2. build per-motif RMSD files (.sh), taking the lowest RMSD of each fragment
3. gather atomtype coordinates (.sh that invokes .py)
4. build distance arrays
(TODO: merge distance arrays of multiple complexes)
then:
5. make histograms
6. (optional)
   score with histograms (single atomtype-atomtype contact)
   OR: 
   score with all histograms (all atomtype-atomtype contacts)

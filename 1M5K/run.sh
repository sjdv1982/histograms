# must already be done: 
#  lrmsd-all.py, build-motif-rmsd-files.sh 
# not used: ../scripts/gather-atomtype-coordinates.sh
# instead, invoke for one case:

mkdir -p coordinates
python3 ../scripts/gather-atomtype-coordinates.py CAC/confr-1.pdb CAC.dat CACr.list coordinates/CAC
mkdir -p distance-arrays
python3 ../scripts/build-distance-arrays.py . --distances 3 10 --discretizations 3 0.5 --rmsd-thresholds 5 15 --rank-chunk 99999999
mkdir -p histograms
python3 ../scripts/make-histograms.py distance-arrays 50 histograms/
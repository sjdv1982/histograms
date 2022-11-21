# - reads ligand Numpy coordinates ($motif-*.npy) from current directory
# - iterates implicitly over all motifs
# - creates cache-$motif/ directory with cache files
set -u -e
currdir=`python3 -c 'import os,sys;print(os.path.dirname(os.path.realpath(sys.argv[1])))' $0`

histo_dir=$1
for i in ???-*.npy; do 
    motif=`echo $i | awk '{print substr($1, 1, 3)}'`
    atomtype=`echo $i | awk '{print substr($1, 5, length($1) - 8)}'`
    echo $motif $atomtype
    mkdir -p cache-$motif
    python3 $currdir/../scripts/discretize_coordinates.py $i $atomtype $histo_dir/1-$atomtype.json cache-$motif
done
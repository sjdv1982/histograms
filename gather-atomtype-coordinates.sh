for c in `cat list`; do
  mkdir $c/coordinates
  for motif in `cat $c/motif.list`; do
    python3 new/gather-atomtype-coordinates.py \
      `head -1 library/$motif.list` \
      $c/$motif-rec.dat \
      library/$motif.list $c/coordinates/$motif
  done
done
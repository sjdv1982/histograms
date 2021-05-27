for c in `cat list`; do
  cd $c
  for motif in `awk '{print $2}' boundfrag.list | sort | uniq`; do
    t1=`mktemp`
    t2=`mktemp`
    t3=`mktemp`
    touch $t1 $t2 $t3
    for frag in `awk -v motif=$motif '$2==motif{print $1}' boundfrag.list`; do
      awk '{print $2}' frag$frag.rmsd > $t2
      paste $t1 $t2 > $t3
      cat $t3 > $t1
    done
    awk '{m=$1; for (n=1;n<=NF;n++) if ($n < m) m=$n; print(m)}' $t1 > $motif.rmsd
    rm -f $t1 $t2 $t3
    echo $c $motif
  done
  cd ..
done
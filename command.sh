awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

bioawk -c fastx '{print $name"\t"length($seq)}' references/Big.MysteryGenomePilon3.PILON.fasta >chrom.sizes

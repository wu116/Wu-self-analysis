#genome prefix
prefix=$1
nearspecies=$2

#generate Allele.ctg.table to identity allele contig
gmap_build -D . -d $prefix genome.nextpolish.fasta
gmap -D . -d $prefix -t 56 -f 2 -n 4 $nearspecies.cds > gmap.gff3

gmap2AlleleTable.pl $nearspecies.gff3

#partition
partition_gmap.pl -g Allele.ctg.table -r genome.nextpolish.fasta -b sample.clean.bam -d wrk_dir

#result: according to the number of chromosome in near species, the draft contig will be partited into several groups.

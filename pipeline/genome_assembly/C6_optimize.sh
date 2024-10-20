prefix=$1

for i in {1..20};
do
	cd wrk_dir/${prefix}${i}
	allhic extract prunning.sub.bam seq.fasta --RE GATC
	allhic optimize group$i.txt prunning.sub.clm

	cd ../../
done

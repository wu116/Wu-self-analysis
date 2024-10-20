prefix=$1

for i in {1..20};
do
	cd wrk_dir/${prefix}${i}

	ALLHiC_rescue -r seq.fasta -b prunning.sub.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt

	cd ../../
done

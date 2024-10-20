prefix=$1

for i in {1..20};
do
	cd wrk_dir/${prefix}${i}

	#for chimeric diploid genome: k = 1
	#for haplotype-resolve genome: k = 2
	ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 1 -m 25

	cd ../../
done

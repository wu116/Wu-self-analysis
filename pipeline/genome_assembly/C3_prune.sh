prefix=$1

for i in {1..20};
do
	cd wrk_dir/${prefix}${i}
	ALLHiC_prune -i ../../Allele.ctg.table -b ../../sample.clean.bam -r seq.fasta
	
	#clean big files
	rm log.txt
	rm removedb_Allele.txt
	rm removedb_nonBest.txt
	rm prunning,sam

	cd ../../
done


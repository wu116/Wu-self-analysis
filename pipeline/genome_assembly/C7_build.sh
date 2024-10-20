prefix=$1

for i in {1..20};
do
	cd wrk_dir/${prefix}${i}
	ALLHiC_build seq.fasta
	cd ../../
done

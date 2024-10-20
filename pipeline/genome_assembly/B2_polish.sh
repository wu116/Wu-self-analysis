fq=$1
prefix=$2

#medaka
medaka_consensus -i $fq -d bridged_contigs.fasta -o ${prefix}_medaka -t 56

#3 rounds racon
minimap2 -t 56 ${prefix}_medaka/consensus.fasta $fq > ${prefix}_racon/round_1.paf 
racon -t 56 $fq ${prefix}_racon/round_1.paf ${prefix}_medaka/consensus.fasta > ${prefix}_racon/racon_round1.fasta

minimap2 -t 56 B13_racon/racon_round1.fasta B13.filter.fastq.gz > B13_racon/round_2.paf 
racon -t 56 B13.filter.fastq.gz B13_racon/round_2.paf B13_racon/racon_round1.fasta> B13_racon/racon_round2.fasta

minimap2 -t 56 B13_racon/racon_round2.fasta B13.filter.fastq.gz > B13_racon/round_3.paf
racon -t 56 B13.filter.fastq.gz B13_racon/round_3.paf B13_racon/racon_round2.fasta> B13_racon/racon_round3.fasta


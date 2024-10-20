#build index of result of nextpolish
bwa index -a bwtsw genome.nextpolish.fasta 
samtools faidx genome.nextpolish.fasta

#map HIC read to result of nextpolish
bwa aln -t 56 genome.nextpolish.fasta clean_combined.1.fq.gz > clean_combined.1.sai 
bwa aln -t 56 genome.nextpolish.fasta clean_combined.2.fq.gz > clean_combined.2.sai
bwa sampe genome.nextpolish.fasta clean_combined.1.sai clean_combined.2.sai clean_combined.1.fq.gz clean_combined.2.fq.gz > sample.bwa_aln.sam

#process and filter the HIC reads
PreprocessSAMs.pl sample.bwa_aln.sam genome.nextpolish.fasta MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam
samtools view -bt genome.nextpolish.fasta.fai sample.clean.sam > sample.clean.bam


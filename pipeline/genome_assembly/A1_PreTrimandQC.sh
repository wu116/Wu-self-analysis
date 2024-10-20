fq=$1

#Pre QC
NanoPlot --fastq $fq -t 20 --N50 -o ${fq}_nanoplot
nanoQC $fq -o ${fq}_nanoQC

#Filter and Trim
pigz -d -p 28 -c $fq|chopper -l 2000 -q 7 --headcrop 50 --tailcrop 20 -t 28|pigz  > $fq.filter.fastq.gz

#Post QC
NanoPlot --fastq $fq.filter.fastq.gz -t 20 --N50 -o ${fq}_filter_nanoplot
nanoQC $fq.filter.fastq.gz -o ${fq}_filter_nanoQC

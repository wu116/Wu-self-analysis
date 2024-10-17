#!/bin/bash

fq=$1
genome=$2
gff=$3
suffix=${4-fq}
length=${5-150}
thread=${6-8}

if [ -z "$fq" ] || [ -z "$genome" ] || [ -z "$gff" ]; then
	echo "Usage: $0 fq genome gff suffix length thread"
	echo "fq: Dir with only sequence raw data files"
	echo "genome: The reference genome fasta file"
	echo "gff: The reference gene annotation gff3 file"
	echo "suffix: The file feature of raw data [fq/fastq]"
	echo "length: The length of the sequence reads [150]"
	echo "thread: CPU [8]"
	exit 1
fi

if ls $fq|grep -vq ".$suffix.gz"; then
	echo "Error: Some files in $fq do not be end with .$suffix.gz"
	exit 1
fi

if ls $genome.*.ht2 1> /dev/null 2>&1;then
	:
else
	hisat2-build -p $thread $genome $genome
fi

[ ! -d  clean ] && mkdir clean
[ ! -d  bam ] && mkdir bam
[ ! -d  quant ] && mkdir quant

for prefix in $(ls -1 $fq|while read line;do echo ${line%[1-2].*};done|sort -u);
do
	~/Music/fastp -w $thread -i ./$fq/${prefix}1.${suffix}.gz \
                -I ./$fq/${prefix}2.${suffix}.gz \
                -o ./clean/${prefix}1_clean.${suffix}.gz \
                -O ./clean/${prefix}2_clean.${suffix}.gz \
                -h ./clean/${prefix}.html -j ./clean/${prefix}.json

	hisat2 --dta -p $thread -x $genome \
                -1 ./clean/${prefix}1_clean.${suffix}.gz \
                -2 ./clean/${prefix}2_clean.${suffix}.gz \
                --summary-file ./bam/${prefix}_hisat2.summary \
                -S ./bam/${prefix}.sam

	samtools sort -@ $thread \
		-o ./bam/${prefix}.sorted.bam \
		./bam/${prefix}.sam

	rm ./bam/${prefix}.sam

	stringtie -e -B -p $thread \
		-G $gff \
		-A ./quant/$prefix/$prefix"_abund".tab \
		-o ./quant/$prefix/$prefix.gtf \
		./bam/${prefix}.sorted.bam
done

prepDE.py -i ./quant -l $length

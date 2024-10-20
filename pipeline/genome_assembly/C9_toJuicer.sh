python path/to/juicebox_scripts/juicebox_scripts/agp2assembly.py groups.agp out.assembly
samtools sort -@ 56 -n sample.bwa_aln.bam > sample.bwa_aln.nSorted.bam

matlock bam2 juicer sample.bwa_aln.nSorted.bam out.links.txt
sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
bash /path/to/3d-dna-201008/visualize/run-assembly-visualizer.sh -p true out.assembly out.sorted.links.txt

rm out.links.txt

#use out.assembly and out.hic to manually correct the assembly by Juicebox

#after correct
#python path/to/juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py \
#	-a ./out.review.ssembly \
#	-f ./genome.nextpolish.fasta \
#	-p review -s

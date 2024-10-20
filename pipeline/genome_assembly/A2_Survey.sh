fq_1=$1
fq_2=$2
prefix=$3

#jellyfish
jellyfish count -C -m 19 -s 20000M -o ${prefix}_km19 -t 56 <(zcat $1) <(zcat $2)
jellyfish histo -t 56 ${prefix}_km19 > ${prefix}_km19.tsv

jellyfish count -C -m 17 -s 20000M -o ${prefix}_km17 -t 56 <(zcat $1) <(zcat $2)
jellyfish histo -t 56 ${prefix}_km17 > ${prefix}_km17.tsv

jellyfish count -C -m 21 -s 20000M -o ${prefix}_km21 -t 56 <(zcat $1) <(zcat $2)
jellyfish histo -t 56 ${prefix}_km21 > ${prefix}_km21.tsv

#genomescope
genomescope.R -i ${prefix}_km19.tsv -o ${prefix}_km19_p2 -k 19 -p 2

genomescope.R -i ${prefix}_km17.tsv -o ${prefix}_km17_p2 -k 17 -p 2

genomescope.R -i ${prefix}__km21.tsv -o ${prefix}_km21_p2 -k 21 -p 2

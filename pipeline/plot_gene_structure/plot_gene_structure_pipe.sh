query=$1
hmm=$2
db=$3
out=$4
evalue=$5
cdd_db=~/database/CDD/Cdd
cdd_model="std"
rpsbproc_data=~/Music/ncbi-blast-2.11.0+/bin/RpsbProc-x64-linux/utils/data
allgff=db/all.gff

#mkdir $out

blastp -query $query -db $db -out $out/$out.blastp.tsv -outfmt 6 -evalue $evalue
hmmsearch --noali -E $evalue --tblout $out/$out.hmm.tsv $hmm $db

cut -f2 $out/$out.blastp.tsv |sort -u > $out/$out.blastp.list
grep -v "^#" $out/$out.hmm.tsv |awk -F ' ' '{print $1}'|sort -u > $out/$out.hmm.list
cat $out/$out.blastp.list $out/$out.hmm.list |sort |uniq -d > $out/$out.final.list

##Man check process

mkdir $out/ManCheck

#CDD
seqkit grep -f $out/$out.final.list $db > $out/$out.final.pep

cat $query $out/$out.final.pep > $out/ManCheck/QueSub.pep

rpsblast -query $out/ManCheck/QueSub.pep -outfmt 11 -evalue 0.01 -out $out/ManCheck/QueSub.cdd.asn -db $cdd_db -sorthsps 2
rpsbproc -i $out/ManCheck/QueSub.cdd.asn -o $out/ManCheck/QueSub.cdd.out -e 0.01 -m $cdd_model -d $rpsbproc_data

#MEME
meme $out/ManCheck/QueSub.pep -o $out/ManCheck/QueSub_meme -nmotifs 10

grep -f<(grep ">" $out/ManCheck/QueSub.pep|sed 's/>//g'|sed 's/ .*//') $allgff > $out/ManCheck/QueSub.gff

python aln2structureplot.py -g $out/ManCheck/QueSub.gff -c $out/ManCheck/QueSub.cdd.out -d 'pfam' -o $out/ManCheck/QueSub.ggtranscript.gff

#align
mafft --auto $out/ManCheck/QueSub.pep > $out/ManCheck/QueSub.aln

#trim
trimal -in $out/ManCheck/QueSub.aln -out $out/ManCheck/QueSub.trim.aln -automated1


##read QueSub.ggtranscript.gff in R to plot gene structure
##use QueSub.trim.aln to plot gene evolution tree

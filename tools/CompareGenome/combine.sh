query=$1
hmm=$2
db=$3
out=$4
evalue=$5

mkdir $out

blastp -query $query -db $db -out $out/$out.blastp.tsv -outfmt 6 -evalue $evalue
hmmsearch --noali -E $evalue --tblout $out/$out.hmm.tsv $hmm $db

cut -f2 $out/$out.blastp.tsv |sort -u > $out/$out.blastp.list
grep -v "^#" $out/$out.hmm.tsv |awk -F ' ' '{print $1}'|sort -u > $out/$out.hmm.list
cat $out/$out.blastp.list $out/$out.hmm.list |sort |uniq -d > $out/$out.final.list


prefix=$1

find wrk_dir/${prefix}* -maxdepth 0|while read line;do
	agp=$line/groups.agp
	grep -v '#' ${agp} | grep 'group' \
		|cut -f 1 |sort -u|while read chrID;do
		ctgNum=`grep -w ${chrID} ${agp} | grep 'W' |wc -l`;
		ctg=`grep -w ${chrID} ${agp} | grep -w  'W' | cut -f 6 | tr '\n' ' '`;
		chrID=`echo ${chrID}|sed 's/group//g'`
	echo -e "${line}_${chrID}\t${ctgNum}\t${ctg}" >> clusters.txt;
	done
done

ALLHiC_rescue -r genome.nextpolish.fasta -b sample.clean.bam -c clusters.txt -i sample.clean.counts_GATC.txt

#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input_tsv_file.tsv>";
die($usage) unless (@ARGV == 1);

my $wrkdir = "wrkdir";
my $input_tsv_file = $ARGV[0];
my @seqname;

system "mkdir $wrkdir";

open(my $in_tsv,"<",$input_tsv_file);
while(<$in_tsv>){
    chomp;
    my @line = split(/\t/);
    my $line_num = $.;
    system "seqkit grep -n -p $line[0] ../../pep_data_mod/O.thomaeum.fa > $wrkdir/$line_num.fa";
    system "seqkit grep -n -p $line[1] ../../../09Synteny/new_1130/subgenome/Zjnewfixed.pep >> $wrkdir/$line_num.fa";
    system "seqkit grep -n -p $line[2] ../../../09Synteny/new_1130/subgenome/Zjnewfixed.pep >> $wrkdir/$line_num.fa";
    system "seqkit grep -n -p $line[3] ../../../09Synteny/new_1130/subgenome/Zjnewfixed.pep >> $wrkdir/$line_num.fa";
    system "seqkit grep -n -p $line[4] ../../../09Synteny/new_1130/subgenome/Zjnewfixed.pep >> $wrkdir/$line_num.fa";
    if(eof $in_tsv) {
        my $Ot = substr($line[0], 0, 7);
        my $Zj_1 = substr($line[1], 0, 4);
        my $Zj_2 = substr($line[2], 0, 4);
        my $Zj_3 = substr($line[3], 0, 4);
        my $Zj_4 = substr($line[4], 0, 4);
        @seqname = ($Ot,$Zj_1,$Zj_2,$Zj_3,$Zj_4);
    }
}
close($in_tsv);

system "ls -1 $wrkdir/*.fa | parallel -j 10 -I {} 'mafft --auto --quiet {} > {.}.aln'";
system "ls -1 $wrkdir/*.aln | parallel -j 10 -I {} '~/Music/trimal-1.4.1/source/trimal -automated1 -in {} -out {.}.trim'";
system "ls -1 $wrkdir/*.trim | parallel -j 10 -I {} 'seqkit seq -w 0 {} > {.}.oneline'";
system "paste -d '' $wrkdir/*.oneline > $wrkdir/all.fa";

open(my $all_fa,"<","$wrkdir/all.fa");
open(my $all_phy,">","$wrkdir/all.phy");
while(<$all_fa>){
    if(m/^>/){
        my $seq = <$all_fa>;
        print $all_phy ">".shift(@seqname)."\n";
        print $all_phy $seq."\n";
    }
}
close($all_phy);
close($all_fa);

system "java -jar ~/Music/prottest-3.4.2/prottest-3.4.2.jar -i $wrkdir/all.phy -all-distributions -F -AIC -BIC -tc 0.5 -threads 10 -o $wrkdir/prottest.out";
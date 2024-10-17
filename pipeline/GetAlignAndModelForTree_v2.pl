#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input_tsv_file.tsv> <fa_dir_with_each_species>";
die($usage) unless (@ARGV == 2);

my $wrkdir = "wrkdir";
my $fa_dir = $ARGV[1];
my $input_tsv_file = $ARGV[0];
my @seqname;
my $sp_num;

if(-d $wrkdir){
    system "rm $wrkdir/*";
} else {
    system "mkdir $wrkdir";
}

open(my $in_tsv,"<",$input_tsv_file);
while(<$in_tsv>){
    chomp;
    my @line = split(/\t/);
    my $line_num = $.;
    if($line_num == 1){
        $sp_num = scalar(@line) - 1;
        @seqname = @line[3 .. $sp_num];
        next;
    }
    for my $k (@seqname){
        print($k)
    }
    for my $i (3 .. $sp_num){
        system "seqkit grep -p $line[$i] $fa_dir/$seqname[$i - 3].fa | seqkit seq -i >> $wrkdir/$line_num.fa";
        die "Error: Seqkit failed.\nseqkit grep -p $line[$i] $fa_dir/$seqname[$i - 3].fa\n" if ($? != 0);
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
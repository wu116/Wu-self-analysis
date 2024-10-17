#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <promoter_pos.txt> <gene_meth.tsv> <threshold>";
die($usage) unless (@ARGV == 3);

#my $promoter_pos_file = "E:/Study/Bioinformatics/genome/temp/newfixed/PROG1/promoter/ZjPROG1_05G020330_promoter_analysis.pos.txt";
#my $gene_meth_file = "E:/Study/Bioinformatics/genome/temp/newfixed/PROG1/promoter/ZjPROG1.tsv2";
#my $threshold = 0.5;
my $promoter_pos_file = $ARGV[0];
my $gene_feature = $1 if $promoter_pos_file =~ m/_(\d+G\d+)_/;
my $gene_meth_file = $ARGV[1];
my $output_file = $promoter_pos_file =~ s/\.pos\.[^.]*$//r;
my $output_pos_file = $output_file.".meth.pos.txt";
my $output_fre_file = $output_file.".meth.fre.txt";
my $threshold = $ARGV[2];

my %fam_pos;
open(my $promoter_pos,"<",$promoter_pos_file) or die "$!";
while(<$promoter_pos>) {
    next if $. == 1;
    chomp;
    my @line = split(/\t/);
    push @{$fam_pos{$line[0]}}, [$line[1], $line[2]];
}
close($promoter_pos);


open(my $gene_meth,"<",$gene_meth_file) or die "$!";
while(<$gene_meth>) {
    next if $. == 1;
    chomp;
    my @line = split(/\t/);
    next if $line[7] < $threshold;
    next unless $line[8] =~ /$gene_feature/;
    my $distance = $line[6];
    foreach my $key (sort keys %fam_pos) {
        foreach my $value (sort {$a->[0] <=> $b->[0]} @{$fam_pos{$key}}) {
            @$value[2] = "non_meth" unless defined @$value[2];
            if ($distance >= @$value[0] && $distance <= @$value[1] && @$value[2] ne "meth"){
                @$value[2] = "meth";
                next;
            }
            elsif (@$value[2] eq "meth") {
                next;
            }
            else {
                @$value[2] = "non_meth";
            }
        }
    }
}
close($gene_meth);

open(my $output_pos,">",$output_pos_file) or die "$!";
open(my $output_fre,">",$output_fre_file) or die "$!";
print $output_pos "fam"."\t"."pos_start"."\t"."pos_end"."\t"."meth"."\n";
print $output_fre "fam"."\t"."total_fre"."\t"."meth_fre"."\n";
foreach my $key (sort keys %fam_pos) {
    my $m = 0;
        foreach my $value (sort {$a->[0] <=> $b->[0]} @{$fam_pos{$key}}) {
            my @out = ($key, @$value[0], @$value[1], @$value[2]);
            print $output_pos join("\t", @out)."\n";
            $m += 1 if @$value[2] eq "meth";
        }
    print $output_fre $key."\t".scalar(@{$fam_pos{$key}})."\t".$m."\n";
}
close($output_fre);
close($output_pos);
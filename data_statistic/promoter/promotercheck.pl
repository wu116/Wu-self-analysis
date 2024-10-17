#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <promoter_PAN.txt>";
die($usage) unless (@ARGV == 1);

my $promoter_PAN_file = $ARGV[0];
#my $promoter_PAN_file = "E:/Study/Bioinformatics/genome/temp/newfixed/PROG1/promoter/ZjPROG1_05G020330_promoter_analysis.txt";
my $output_file = $promoter_PAN_file =~ s/\.[^.]*$//r;
my $output_pos_file = $output_file.".pos.txt";
my $output_fre_file = $output_file.".fre.txt";

my %fam_id;

open(my $promoter_PAN,"<",$promoter_PAN_file) or die "$!";
while(<$promoter_PAN>) {
    chomp;
    next if $. == 1;
    my @line = split(/\t/);
    next if $line[1] =~ /\(.*\)/;
    my @fam = split(/;/, $line[1]);
    my $pos = $line[2];
    my $seq = $line[5];
    foreach my $i (@fam) {
        $i =~ s/(^\s+|\s+$)//g;
        if (not exists $fam_id{$i}) {
            $fam_id{$i} = [[$pos, length($seq)]];
        }
        elsif (exists $fam_id{$i}){
            my @new_fam_id = @{$fam_id{$i}};
            my $n = 0;
            foreach my $j (@new_fam_id) {
                if (abs($pos - @$j[0]) <= 5){
                    last;
                }
                else {
                    $n += 1;
                }
            }
            if ($n == scalar(@new_fam_id)) {
                push @{$fam_id{$i}}, [$pos, length($seq)];
            }
        }
    }
}
close($promoter_PAN);

open(my $output_pos,">",$output_pos_file);
print $output_pos "fam"."\t"."pos_start"."\t"."pos_end"."\n";
foreach my $key (sort keys %fam_id) {
    foreach my $value (sort {$a->[0] <=> $b->[0]} @{$fam_id{$key}}) {
        my @out = ($key, @$value[0], @$value[0] + @$value[1] - 1);
        print $output_pos join("\t", @out)."\n";
    }
}
close($output_pos);

open(my $output_fre,">",$output_fre_file);
print $output_fre "fam"."\t"."fre"."\n";
foreach my $key (sort keys %fam_id) {
    my @out = ($key, scalar(@{$fam_id{$key}}));
    print $output_fre join("\t", @out)."\n";
}
close($output_fre);
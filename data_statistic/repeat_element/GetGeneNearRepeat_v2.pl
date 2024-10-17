#!/usr/bin/perl
use strict; use warnings; use Getopt::Std;

my $usage = "$0 -f <check_feature(gene/mRNA)> \
    -l <left_bp> -r <right_bp> \
    -G <check_file.gff3> -R <repeat_file.gff3> -O <output_file.gff3>";

my %opt;
getopts('f:l:r:G:R:O:', \%opt);
my ($chr,$feature,$attributes,$start,$end,$strand);
my $match;

open(OUTPUT,">","$opt{O}");
print "";
close(OUTPUT);

open(CF,"$opt{G}") or die $!;
open(my $RF,"$opt{R}") or die $!;
while(<CF>){
    next if (m/^#/);
    next if (m/^$/);
    ($chr,undef,$feature,$start,$end,undef,$strand,undef,$attributes) = split(/\t/);
    next unless ($feature eq $opt{f});
    $match = 0;
    my $fname = $1 if($attributes =~ /ID=(.*?)(?:;|$)/);
    my $clbp = $opt{l};
    my $crbp = ($start - $end + 1 < $opt{r} ? $start - $end + 1 : $opt{r});
    my $RF_FH = 0;

    open(OUTPUT,">>","$opt{O}");
    seek($RF,$RF_FH,0);
    while(<$RF>){
        next if (m/^#/);
        next if (m/^$/);
        my ($RF_chr,undef,$RF_feature,$RF_start,$RF_end,undef,$RF_strand,undef,$RF_attributes) = split(/\t/);
        next unless ($chr eq $RF_chr);
        if($start > $RF_start + 2 * $clbp - 1){
            $RF_FH = tell $RF;
            next;
        }

        if($strand eq "+"){
            if($start <= $RF_end + $clbp - 1 and $start >= $RF_start - $crbp + 1){
                print OUTPUT join("\t",($chr,$start,$end,$strand,$fname,$RF_feature,$RF_start,$RF_end,$RF_strand,$RF_attributes));
                $match = 1;
            }
        }elsif($strand eq "-"){
            if($end <= $RF_end + $crbp - 1 and $end >= $RF_start - $clbp + 1){
                print OUTPUT join("\t",($chr,$start,$end,$strand,$fname,$RF_feature,$RF_start,$RF_end,$RF_strand,$RF_attributes));
                $match = 1;
            }
        }

         if ($match == 1 and $end < $RF_start - 2 * $clbp + 1){
             last;
         }
    }
    close(OUTPUT);
}
close($RF);
close(CF);


sub merge_sort_recursive{
    my $array = @_;

    return $array if @$array <= 1;

    my $middle = int(@$array);
    my @left = merge_sort_recursive([@$array[0..($middle-1)]]);
    my @right = merge_sort_recursive([@$array[$middle..$#$array]]);
    return merge(\@left,\@right);
}
sub merge{
    my ($left,$right) = @_;

    my @result;
    while(@$left && @$right){
        my $left_start = (split(/\t/,$left->[0]))[3];
        my $left_end = (split(/\t/,$left->[0]))[4];
        my $right_start = (split(/\t/,$right->[0]))[3];
        my $right_end = (split(/\t/,$right->[0]))[4];

        if($left_start < $right_start || ($left_start == $right_start && $left_end <= $right_end)){
            push(@result,shift(@$left));
        }else{
            push(@result,shift(@$right));
        }
    }
    push(@result,@$left,@$right);
    return(@result);
}
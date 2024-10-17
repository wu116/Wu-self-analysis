#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my $usage = "$0 -c <chr.agp> \
    -A <Allele.ctg.table> \
    -o <redundantAllelecontig.tsv>";

my %opt;
getopts('c:A:o:', \%opt);

my %anchored_data;

open(my $chr_agp, "<", "$opt{c}") or die "Can not open input file $opt{c}: $!";
while (my $line = <$chr_agp>){
    next if $line =~ /^#/;
    chomp $line;
    my @split_line = split(/\t/, $line);
    my ($anchored_chr_id, $type, $anchored_contig_id) = @split_line[0, 4, 5];
    next if $type eq "U";
    $anchored_chr_id =~ s/chr/Chromosome/;

    $anchored_data{$anchored_contig_id} = $anchored_chr_id;
}
close($chr_agp);

my %allele_data;

open(my $Allele_ctg_table, "<", "$opt{A}") or die "Can not open input file $opt{A}: $!";
while (my $line = <$Allele_ctg_table>){
    chomp $line;
    my @split_line = split(/\t/, $line);
    my $chr_id = $split_line[0];
    my @ctg_ids = @split_line[2 .. $#split_line];
    next if @ctg_ids < 2;
    foreach my $ctg_id (@ctg_ids){
        if (exists($anchored_data{$ctg_id}) && $anchored_data{$ctg_id} eq $chr_id) {
            my @allele_ctg_ids = grep {$_ ne $ctg_id && !exists($anchored_data{$_})} @ctg_ids;
            if (@allele_ctg_ids > 0){
                $allele_data{$ctg_id} = [] unless exists($allele_data{$ctg_id});
                foreach my $allele_ctg_id (@allele_ctg_ids){
                    push(@{$allele_data{$ctg_id}}, $allele_ctg_id) unless grep {$allele_ctg_id eq $_} @{$allele_data{$ctg_id}};
                }
            }
        }
    }
}
close($Allele_ctg_table);

open(my $output, ">", "$opt{o}") or die "Can not open output file $opt{o}: $!";
foreach my $key (keys (%allele_data)){
    print $output $key."\t";
    my $array = $allele_data{$key};
    foreach my $element (@$array){
        print $output $element."\t";
    }
    print $output "\n";
}
close($output);

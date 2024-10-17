#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input.pgenes.align.gz> <output.gff3>";
die ($usage) unless (@ARGV == 2);

my $input_align = $ARGV[0];
my $output_gff3 = $ARGV[1];
my $source = "PseudoPipe";

open (my $input, "<", $input_align) or die "Can not open input file $input_align: $!";
open (my $output, ">", $output_gff3) or die "Can not open output file $output_gff3: $!";
while (my $line = <$input>){
    next unless $line =~ /^>/;
    chomp $line;
    my @split_line = split (/\s{2}/, $line);
    my ($pseudo_name, $start, $end, $strand, $query, $exon_bound) = @split_line[0, 2, 3, 4, 5, 19];
    my $pseudo_id = $1."===$start" if $pseudo_name =~ />(.*)/;
    my $chr_id = $1 if $pseudo_name =~ />(.*?)===/;

    my @output_pseudogene = ($chr_id, $source, "pseudogene", $start, $end, ".", $strand, ".", "ID=$pseudo_id");
    print $output join("\t", @output_pseudogene)."\n";

    if ($exon_bound =~ /^com/){
        my @exon_group = split (/\s/, $exon_bound);
        my $i = 0;
        foreach my $exon_element (@exon_group){
            my ($exon_start, $exon_end) = ($1, $2) if $exon_element =~ /(\d+)\.\.(\d+)/;
            my @output_exon = ($chr_id, $source, "exon", $exon_start, $exon_end, ".", "-", ".", "ID=$pseudo_id.exon$i;Parent=$pseudo_id");
            print $output join("\t", @output_exon)."\n";
        }
    }elsif ($exon_bound =~ /^\(/){
        my @exon_group = split (/\s/, $exon_bound);
        my $i = 0;
        foreach my $exon_element (@exon_group){
            my ($exon_start, $exon_end) = ($1, $2) if $exon_element =~ /(\d+)\.\.(\d+)/;
            my @output_exon = ($chr_id, $source, "exon", $exon_start, $exon_end, ".", "+", ".", "ID=$pseudo_id.exon$i;Parent=$pseudo_id");
            print $output join("\t", @output_exon)."\n";
        }
    }
}
close ($output);
close ($input);
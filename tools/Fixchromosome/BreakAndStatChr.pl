#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input_chr.fa> <output_chr.fa> <broken_contig.list>";
die($usage) unless (@ARGV == 2);

my $fasta_input_file = $ARGV[0];
my $broken_contig_file = $ARGV[1];

my $seq_id;
my $seq = "";
my @seq_data;

open(my $fasta_input,"<",$fasta_input_file) or die "$!";
while(<$fasta_input>){
    chomp;
    $seq_id = $1 and next if (m/^>(.*)/);
    $seq .= $_;
}
close($fasta_input);

@seq_data = $seq =~ /([ACGT]+|N+)/g;
my $num = 1;
my $length;
my $start = 0;
my $end = 0;

open(my $broken_contig,">",$broken_contig_file);
for my $sub_seq (@seq_data){
    $start = $end + 1;
    $length = length($sub_seq);
    $end = $start + $length - 1;
    if ($sub_seq =~ /[ACGT]+/){
        my @output = ($seq_id,"contig".$num,$start,$end,$length,"W");
        print $broken_contig join("\t",@output)."\n";
        $num += 1;
    }elsif ($sub_seq =~ /N+/){
        my @output = ($seq_id,"contig".$num,$start,$end,$length,"U");
        print $broken_contig join("\t",@output)."\n";
        $num += 1;
    }else{
        die "Split seq failed: $!";
    }
}
close($broken_contig);
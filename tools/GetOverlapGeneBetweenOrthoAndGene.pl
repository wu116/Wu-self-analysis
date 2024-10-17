#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <input.ortho.txt> <input.gene.list> <output.ortho.txt>";
die($usage) unless (@ARGV == 3);

my $input_ortho_file = $ARGV[0];
my $input_gene_list = $ARGV[1];
my $output_ortho_file = $ARGV[2];

my %ortho_hash;
my @gene_list;

open(my $in_list,"<",$input_gene_list);
while(<$in_list>){
    chomp;
    push(@gene_list,$_);
}
close($in_list);

open(my $in_ortho,"<",$input_ortho_file);
while(my $line = <$in_ortho>){
    if($line =~ /^([^:]+): (.+)$/){
        chomp $line;
        my $orthogroup_id = $1;
        my @gene_id_array = split(/\s+/,$2);
        $ortho_hash{$orthogroup_id} = \@gene_id_array;
    }
}
close($in_ortho);

open(my $out_ortho,">",$output_ortho_file);
my %filter_ortho_hash = FilterArrayInHash(\%ortho_hash,\@gene_list);
while(my ($key,$value) = each(%filter_ortho_hash)){
    print $out_ortho $key.": ".join(" ",@$value)."\n";
}
close($out_ortho);

sub FilterArrayInHash{
    my ($hash_sub,$array_sub) = @_;
    my %hash = %$hash_sub;
    my @ref_array = @$array_sub;
    my %filter_hash;
    foreach my $key (keys %hash){
            my @value = @{$hash{$key}};
            #my @filter_value = grep {$_ ~~ @ref_array} @value;
            my @filter_value;
            foreach my $val (@ref_array){
                if(grep {$_ eq $val} @value){
                    push(@filter_value,$val);
                }
            }
        if(scalar(@filter_value) > 0) {
            $filter_hash{$key} = \@filter_value;
        }
    }
    return(%filter_hash)
}
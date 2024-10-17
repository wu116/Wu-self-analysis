#!/usr/bin/perl


use strict; use warnings;

my $usage = "trf_2_fasta.pl <trf dat file>";
die "$usage" unless (@ARGV == 1);
my $seq_counter = 0;
my $header;

while(<>){
	#skip blank lines
	next if (m/^$/);

	#extract Seq header
	if (m/^Sequence: (.*)/){
		$header = $1;
		$seq_counter++;
	}

	next unless (m/^\d+ \d+ \d+ \d+\.\d /);

	#capture
	my ($start,$end,$period,$copies,$length,$identity,$indels,$score,$a,$c,$g,$t,$entropy,$seq) = split(/\s+/);

	print $header."\t".$start."\t".$end."\t".$header."_".$copies."_".$length."\n";
}

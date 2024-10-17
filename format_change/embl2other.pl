#!/usr/bin/perl
use strict; use warnings; use Getopt::Std;

my $usage = "embl2other.pl -e <embl file> -g <output gff3 file> or -f <output fasta file>";

my %opt;
getopts('e:g:f:', \%opt);

my $AC;

if($opt{f}){
    my $SQ;

	open(FASTA,">","$opt{f}");
	print "";
	close(FASTA);

	open(FASTA,">>","$opt{f}");
	open(EMBL,"$opt{e}") or die $!;
	while(<EMBL>){
	    if(m/^AC/){
		    next if (split(/\s+/))[1] eq "XXX;";
		    $AC = (split(/\s+/))[2];
		    print FASTA ">".$AC."\n";
		}

		if(m/^SQ/ .. m/^\/\//){
			next if (m/^SQ/);
			next if (m/^\/\//);
			s/[^A-Z]//g;
			$SQ = $_;
			print FASTA $SQ."\n";
		}
	}
	close(FASTA);
	close(EMBL);
}

if($opt{g}){
    my $feature;
    my @start;
    my @end;
    my $strand;
    my $LOC;
    my @ID;
    my $source;
    my $gene_id;
    my @mRNA_id;
	my $cds_num;

	open(GFF,">","$opt{g}");
        print "";
        close(GFF);

	open(GFF,">>","$opt{g}");
	open(EMBL,"$opt{e}") or die $!;
	while(<EMBL>){
		if(m/^ID/ .. m/^\/\//){
			if(m/^AC/){
			    next if (split(/\s+/))[1] eq "XXX;";
			    $AC = (split(/\s+/))[2];
			}

			if(m/^FT/){
				next if (m/ source|\/mol_type|\/organism|gap|\/estimated_length/);
				my @FT_line = split(/\s+/);
				if(@FT_line == 3){
					$feature = $FT_line[1];
					my ($start_res,$end_res) = getcoor($_);
					@start = @$start_res;
					@end = @$end_res;
					$strand = ($FT_line[2] =~ /complement/ ? "-" : "+");
					@ID = ();
				}

				if(@FT_line == 2){
					next if ($feature =~ /3'UTR|5'UTR/);
					next if (m/\/transl_table=/);
					$LOC = $1 if(m/\/locus_tag="(.*)"/);
					push(@ID,$1) if(m/\/note="ID:(.*)"/);
					$source = $1 if(m/\/note="source:(\w+)"/);
				}

				if($feature eq "gene" and $source){
					my @gene_output = ($AC,$source,$feature,$start[0],$end[0],"\.",$strand,"\.","ID=$ID[0];locus_tag=$LOC");
					print GFF join("\t",@gene_output)."\n";
					$source = undef;
					$gene_id = $ID[0];
					@mRNA_id = ();
					$cds_num = 0;
				}

				if($feature eq "mRNA" and $source){
					my @mRNA_output = ($AC,$source,$feature,$start[0],$end[-1],"\.",$strand,"\.","ID=$ID[0];Parent=$gene_id");
					print GFF join("\t",@mRNA_output)."\n";
					push(@mRNA_id,$ID[0]);

					for my $e_num (0 .. @start-1){
						my @exon_output = ($AC,$source,"exon",$start[$e_num],$end[$e_num],"\.",$strand,"\.","ID=$ID[0].exon.".($e_num+1).";Parent=$ID[0]");
						print GFF join("\t",@exon_output)."\n";
					}

					$source = undef;
				}

				if($feature eq "exon" and $source){
					$source = undef;
					next;
				}

				if($feature eq "CDS" and $source){
					my @phase = (0);
					if($strand eq "+" and @start > 1){
						for my $num (0 .. @start-2){
							$phase[$num+1] = (3-($end[$num]-$start[$num]+1+$phase[$num])%3);
							$phase[$num+1] = 0 if ($phase[$num+1] == 3);
						}
						for my $i (0 .. @start-1){
							my @CDS_output = ($AC,$source,$feature,$start[$i],$end[$i],"\.",$strand,$phase[$i],"ID=$ID[0].cds.".($i+1).";Parent=$mRNA_id[$cds_num]");
							print GFF join("\t",@CDS_output)."\n";
						}
					}elsif($strand eq "-" and @start > 1){
						for my $num (0 .. @start-2){
							$phase[$num+1] = (3-($end[@start-2-$num+1]-$start[@start-2-$num+1]+1+$phase[$num])%3);
							$phase[$num+1] = 0 if ($phase[$num+1] == 3);
						}
						for my $i (0 .. @start-1){
							my @CDS_output = ($AC,$source,$feature,$start[$i],$end[$i],"\.",$strand,$phase[@start-1-$i],"ID=$ID[0].cds.".($i+1).";Parent=$mRNA_id[$cds_num]");
							print GFF join("\t",@CDS_output)."\n";
						}
					}elsif(@start == 1){
						my @CDS_output = ($AC,$source,$feature,$start[0],$end[0],"\.",$strand,"0","ID=$ID[0].cds."."1".";Parent=$mRNA_id[$cds_num]");
                                                print GFF join("\t",@CDS_output)."\n";
					}
					$source = undef;
					$cds_num += 1;
					die("Number of CDS and Number of mRNA could not be paired!") if ($cds_num > scalar(@mRNA_id));
				}
			}
		}
	}
	close(EMBL);
	close(GFF);
}


sub getcoor{
	my ($line) = @_;
	my @start;
	my @end;
	chomp $line;
	if($line =~ /,$/){
		print $line."\n";
		chomp $line;
		print $line."\n";
		$line = (split(/\s+/,$line))[2];
		my $next_line = <EMBL>;
		while($next_line =~ /,$/){
			chomp $next_line;
			my $next_line_split = (split(/\s+/,$next_line))[1];
			$line .= $next_line_split;
			$next_line = <EMBL>;
		}
		chomp $next_line;
		my $next_line_split = (split(/\s+/,$next_line))[1];
		$line .= $next_line_split;
	}
	my @group = split(",",$line);
	foreach my $group (@group){
		if($group =~ /(\d+)\.\.(\d+)/){
			push(@start,$1);
			push(@end,$2);
		}
	}
	return(\@start,\@end);
}
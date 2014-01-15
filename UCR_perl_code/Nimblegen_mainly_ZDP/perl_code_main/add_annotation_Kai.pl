#!/usr/bin/perl -w
# add annotation to the identified methylated regions
use strict;
use Bio::SeqIO;

my $debug = 0;
my $usage = "$0 <MAT output> <gene_TE GFF> <func annot> <TE frags> <intergenic>";
die $usage unless(@ARGV >= 5);
my ($mat_bed, $gene_TE, $func, $TE_frag, $inter_file) = @ARGV[0..4];
open(MT, $mat_bed) or die "Can't open $mat_bed:$!";
my %mat_regions; # chr->start=\@
while(<MT>){
	next if(/^browser/);
	next if(/^track/);
	next unless(/\w+/);
	chomp;
	my ($chr, $start, $end, $p, $mstart, $mend) = split /\t/;
	if($chr eq "chloroplast"){
		$chr = "chrc";
	}
	if($chr eq "mitochondria"){
		$chr = "chrm";
	}
	$mat_regions{$chr}->{$start} = [$start, $end, $p, $mstart, $mend];
}

my %func;
open(FF, $func) or die "Can't open $func: $!";
while(<FF>){
	next if(/^Model_name/);
	chomp;
	my @temp = split /\t/;
	if($temp[0] =~ /(AT\wG\d+)\.\d+/){
		my $gene = $1;
		if(!defined $func{$gene}){
			if(defined $temp[2] && $temp[2] =~ /\w+/){
			   $func{$gene} = $temp[2];
			}elsif(defined $temp[3] && $temp[3] =~ /\w+/){
				$func{$gene} = $temp[3];
			}elsif(defined $temp[4] && $temp[4] =~ /\w+/){
				$func{$gene} = $temp[4];
			}else{
				$func{$gene} = "NONE";
			}
		}
	}else{
		die "gene name ", $temp[0], " does not match pattern";
	}
}
	

my %genes;
my %notes;
my %annot;
open(GFF, $gene_TE) or die "Can't open $gene_TE:$!";
while(<GFF>){
	chomp;
	my @temp = split /\t/;
	my $chr = lc $temp[0];
	if($temp[2] eq "gene" || $temp[2] eq "pseudogene" || $temp[2] eq "transposable_element_gene"){
		foreach my $st(keys %{$mat_regions{$chr}}){
			my ($start, $end, $p,$mstart, $mend) = @{$mat_regions{$chr}->{$st}};
			if(!($end < $temp[3] || $start > $temp[4])){
				my @tag_vals = split /;/, $temp[8];
				my %vals;
				foreach my $tag_val(@tag_vals){
					my @t = split /=/, $tag_val;
					$vals{$t[0]} = $t[1];
				}
				my $id = "NONE";
				if(defined $vals{"ID"}){
					$id = $vals{"ID"};
				}
				my $note = "NONE";
				if(defined $vals{"Note"}){
					$note = $vals{"Note"};
				}
                
				my $fun = "NONE";
				if(defined $func{$id}){
					$fun = $func{$id};
				}
				if(defined $genes{$chr}->{$start}){
					$genes{$chr}->{$start} .= ";$id";
					$notes{$chr}->{$start} .= ";$note";
					$annot{$chr}->{$start} .= ";$fun";
				}else{
					$genes{$chr}->{$start} = $id;
					$notes{$chr}->{$start} = $note;
					$annot{$chr}->{$start} = $fun;
				}
			}
		}
	}
}

my %TE;
my %TE_family;
open(TF, $TE_frag) or die "Can't open $TE_frag:$!";
while(<TF>){
	next if(/^Transposon_Name/);
	chomp;
	my @temp = split /\t/;
	my $chr = "chr1";
	if($temp[0] =~ /AT(\d)TE/){
		$chr = "chr" . $1;
	}else{
		die "Transposon ", $temp[0], " doesn't match pattern";
	}
	foreach my $st(keys %{$mat_regions{$chr}}){
			my ($start, $end, $p,$mstart, $mend) = @{$mat_regions{$chr}->{$st}};
			if(!($end < $temp[2] || $start > $temp[3])){
				if(defined $TE{$chr}->{$start}){
					$TE{$chr}->{$start} .= ";" . $temp[0];
					$TE_family{$chr}->{$start} .= ";" . $temp[4] . ":" . $temp[5];
				}else{
					$TE{$chr}->{$start} = $temp[0];
					$TE_family{$chr}->{$start} = $temp[4] . ":" . $temp[5];
				}
			}
		}
}

## now intergenic
my $seqin = Bio::SeqIO->new(-file=>$inter_file, -format=>'fasta');
my %intergenic;
while(my $seq = $seqin->next_seq){
	my $chr = "chr1";
	my ($s, $e) = (0,0);
	if($seq->desc =~ /(chr\w):(\d+)\-(\d+)/){
		$chr = lc  $1;
	    ($s, $e) = ($2, $3);
		if($debug){
		    print STDERR "Intergenic start=$s, end=$e\n";
		}
	}else{
		die "Intergenic seq ", $seq->id, " does not have expected description: ", $seq->desc;
	}
	foreach my $st(keys %{$mat_regions{$chr}}){
			my ($start, $end, $p,$mstart, $mend) = @{$mat_regions{$chr}->{$st}};
			if(!($end < $s || $start > $e)){
				if(defined $intergenic{$chr}->{$start}){
					$intergenic{$chr}->{$start} .= ";" . $seq->id;
				}else{
					if($debug){
						print STDERR "Adding intergenic ", $seq->id, "\n";
					}
					$intergenic{$chr}->{$start} = $seq->id;
				}
			}
		}
}
print join("\t", ("Chr", "Peak_Start", "Peak_End", "P","maxfour_start","maxfour_end" ,"Gene", "GeneType", "GeneAnnot",
                  "TE", "TEFamily", "Intergenic")), "\n";
foreach my $chr(sort keys %mat_regions){
	my %reg = %{$mat_regions{$chr}};
    foreach my $start(sort keys %reg){
		my ($start, $end, $p,$mstart, $mend) = @{$reg{$start}};
		my ($gene, $geneType, $geneAnnot, $TE, $TEFamily, $inter) = 
		("NONE", "NONE", "NONE", "NONE", "NONE", "NONE");
		if(defined $genes{$chr}->{$start}){
			$gene = $genes{$chr}->{$start};
		}
		if(defined $notes{$chr}->{$start}){
			$geneType = $notes{$chr}->{$start};
		}
		if(defined $annot{$chr}->{$start}){
			$geneAnnot = $annot{$chr}->{$start};
		}
		if(defined $TE{$chr}->{$start}){
			$TE = $TE{$chr}->{$start};
		}
		if(defined $TE_family{$chr}->{$start}){
			$TEFamily = $TE_family{$chr}->{$start};
		}
		if(defined $intergenic{$chr}->{$start}){
			$inter = $intergenic{$chr}->{$start};
		}
		print join("\t", ($chr, $start, $end, $p,$mstart, $mend, $gene, $geneType, $geneAnnot,
		     $TE, $TEFamily, $inter)), "\n";
	}
}			

print STDERR "\a";
exit;
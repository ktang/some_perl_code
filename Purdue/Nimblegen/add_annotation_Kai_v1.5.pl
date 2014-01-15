#!/usr/bin/perl -w
# add annotation to the identified methylated regions

=head2
version 1.5(from v3)

add strand information
perl add_annotation.pl at1g54840.bed ~/data/ath/tair9/TAIR9_GFF3_genes_transposons.gff ~/data/ath/tair9/TAIR9_functional_descriptions ~/data/ath/tair9/TAIR9_Transposable_Elements.txt ~/data/ath/tair9/TAIR9_intergenic_20090619 >at1g54840_bed_with_annot.txt

use bed input, no trimmean.
=cut

#BEGIN{unshift @INC, "/Library/Perl/5.10.0/BioPerl";}
use strict;
use Bio::SeqIO;

my $debug = 0;
my $usage = "$0 <bed_file>";
die $usage unless(@ARGV == 1);
my $bed= $ARGV[0];



open(MT, $bed) or die "Can't open $bed:$!";

my $func = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_functional_descriptions";
open(FF, $func) or die "Can't open $func: $!";

my $gene_TE = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_GFF3_genes_transposons.gff";
open(GFF, $gene_TE) or die "Can't open $gene_TE:$!";

my $TE_frag = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_Transposable_Elements.txt";
open(TF, $TE_frag) or die "Can't open $TE_frag:$!";

my $inter_file = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_intergenic_20090619";

open(INTER , $inter_file) or die "cannot open $inter_file:$!";
close(INTER);

my %mat_regions; # chr->start=\@
while(<MT>){
	next if(/^browser/);
	next if(/^track/);
	next unless(/\w+/);
	chomp;
	my ($chr, $start, $end) = split /\t/;
	if($chr eq "chloroplast"){
		$chr = "chrc";
	}
	if($chr eq "mitochondria"){
		$chr = "chrm";
	}
	$mat_regions{$chr}->{$start} = [$start, $end];
}

my %func;
#open(FF, $func) or die "Can't open $func: $!";
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
	
######
my %strands;
#####


my %genes;
my %notes;
my %annot;
#open(GFF, $gene_TE) or die "Can't open $gene_TE:$!";
while(<GFF>){
	chomp;
	my @temp = split /\t/;
	my $chr = lc $temp[0];
	if($temp[2] eq "gene" || $temp[2] eq "pseudogene" || $temp[2] eq "transposable_element_gene"){
		foreach my $st(keys %{$mat_regions{$chr}}){
		#	my ($start, $end, $p,$mstart, $mend) = @{$mat_regions{$chr}->{$st}};
			my ($start, $end) = @{$mat_regions{$chr}->{$st}};
		
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
				
				############
				$strands{$chr}->{$start} = $temp[6];
				############
				
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
#open(TF, $TE_frag) or die "Can't open $TE_frag:$!";
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
print join("\t", ("Chr", "Peak_Start", "Peak_End","Gene", "Strand","GeneType", "GeneAnnot",
                  "TE", "TEFamily", "Intergenic")), "\n";
foreach my $chr(sort keys %mat_regions){
	my %reg = %{$mat_regions{$chr}};
    foreach my $start(sort keys %reg){
	#	my ($start, $end, $p,$mstart, $mend) = @{$reg{$start}};
		my ($start, $end) = @{$reg{$start}};
	
		my ($gene, $geneType, $geneAnnot, $TE, $TEFamily, $inter , $strand) = 
		("NONE", "NONE", "NONE", "NONE", "NONE", "NONE" , "NONE");
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
		if(defined $strands{$chr}->{$start}){
			$strand = $strands{$chr}->{$start};
		}
		print join("\t", ($chr, $start, $end, $gene, $strand, $geneType, $geneAnnot,
		     $TE, $TEFamily, $inter)), "\n";
	}
}			


exit;
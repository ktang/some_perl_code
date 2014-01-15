#!/usr/bin/perl -w
=head2
 
 
 not finished!!!!!!!!!
 
 
Dec 1,2010
Kai TANG

with previous annotation, know which gene the 
peak loacate.

read every line and if a gene, check
wheather it exists.

then put all the line relate to this gene 
into array.

CDS
chromosome
exon
five_prime_UTR
gene  # this
mRNA
miRNA
ncRNA
protein
pseudogene  # this
pseudogenic_exon
pseudogenic_transcript
rRNA
snRNA
snoRNA
tRNA
three_prime_UTR
transposable_element_gene #and this

=cut

use strict;
my $debug = 0;
my $usage = "$0 <input_annotation_file> <gff_file>";

die $usage unless (@ARGV == 2);

my ($inFile,$gffFile) = @ARGV[0..1];

open (IN,$inFile)
	or die "cannot open $inFile:$!";
my %peaks_h;#all line;
my %ends_h;
my %genes_h;
my %arr_h;# key = gene_name; value is a array contains all the lines assotiated the gene.

my %intron_h;
my %cds_h;
my %three_h;
my %five_h;


while (<IN>){
	chomp;
	my @pts = split /\t/;
	$peaks_h{$pts[0]}->{$pts[1]}= \@pts;
	$ends_h{$pts[0]}->{$pts[1]} = $pts[2];
	
#	$intron_h{$pts[0]}->{$pts[1]} = 0;
#	$cds_h{$pts[0]}->{$pts[1]} = 0;
#	$three_h{$pts[0]}->{$pts[1]} = 0;
#	$five_h{$pts[0]}->{$pts[1]} = 0;

	if ($pts[4] =~ /;/){
		my ($gene1,$gene2) = split /;/,$pts[4];
		$genes_h{$gene1} = [$pts[0],$pts[1],$pts[2]];
		$genes_h{$gene2} = [$pts[0],$pts[1],$pts[2]];
	}
	else{
		$genes_h{$pts[4]} = [$pts[0],$pts[1],$pts[2]];
	}
}

close(IN);

my %genes_cor_h;#sotre the start and end of all genes

open(GFF,$gffFile)
	or die "cannot open $gffFile:$!";

while (<GFF>){
	next unless (/gene/);
	chomp;
	my @pts = split /\t/;
	my $key;
	if ($pts[8] =~ /^ID=(AT.G.....);/){
		$key = $1;
	}	
	$genes_cor_h{$key} = [$pts[3],$pts[4]];
	my $line;
IF:	if (exists $genes_h{$key}){
		while ($line = <GFF>){
			push @{$arr_h{$key}},$line  unless ($line =~ /gene/);
			last if ($line =~ /gene/);
		}
		
		my @temp = split /\t/,$line;
		if ($temp[8] =~ /^ID=(AT.G.....);/){
			$key = $1;
			goto IF;
		}

	}		
}
my $i = 0;

foreach my $key (sort keys %genes_h){
	$i ++;
	
	my @lines = @{$arr_h{$key}};
	my ($chr,$start,$end) = @{$genes_h{$key}};
	my ($g_start,$g_end) = @{$genes_cor_h{$key}};
	
	if ($debug){
		if ($i <= 3){
			print STDERR "$i:\n",@lines,"\n";
		}
	}
	
	for (my $j = 0; $j <= $#lines;$j++){
			my $this = $lines[$j];
			chomp $this;
			my @pts = split /\t/,$this;
			
			if (overlap ($start,$end, $pts[3],$pts[4])){
				if ($pts[2] eq "CDS"){
					
					if (defined $cds_h{$chr}->{$start}){
					
						my $ori = $cds_h{$chr}->{$start};
						my $post=";".$pts[0].":".$pts[3]."-".$pts[4];
						$cds_h{$chr}->{$start} =  $ori.$post ;
					}
					else{
						$cds_h{$chr}->{$start} = $pts[0].":".$pts[3]."-".$pts[4];
					}
				}
				elsif ($pts[2] eq "five_prime_UTR"){
					if (defined $five_h{$chr}->{$start}){
						my $ori = $five_h{$chr}->{$start};
						my $post = ";".$pts[0].":".$pts[3]."-".$pts[4];
						$five_h{$chr}->{$start} = $ori.$post;
					}
					else{
						
						$five_h{$chr}->{$start} = $pts[0].":".$pts[3]."-".$pts[4];
					}
				}
				elsif ($pts[2] eq "three_prime_UTR"){
					if (defined $three_h{$chr}->{$start} ){
						my $ori = $three_h{$chr}->{$start} ;
						my $post = ";".$pts[0].":".$pts[3]."-".$pts[4];
						$three_h{$chr}->{$start} = $ori.$post;
					}
					else{
						$three_h{$chr}->{$start} = $pts[0].":".$pts[3]."-".$pts[4];
					}
				}

			}
	}
	
	
	my $first_line = $lines[0];
	chomp $first_line;
	my @first_pts = split /\t/,$first_line;
	
	my ($last_cor, $curr_cor ) = (0,0);

	if ($first_pts[6] eq '+'){
		
		my $mRNA_flag = 0;
		for (my $j = 0; $j <= $#lines;$j++){
			
			
			my $this = $lines[$j];
#			next unless ($this =~ /exon|CDS/);
			chomp $this;
			my @pts = split /\t/,$this;
			my ($f_start, $f_end) = @pts[3..4];
			my $m_start = 0;
			
			if ($pts[2] eq 'mRNA'){
				if ($mRNA_flag == 0){
					$mRNA_flag = 1;
				}
				else{
					
				}
					
				$m_start = $f_start;
				if ($f_start > $g_start){
					print STDERR $lines[$j];
					if (overlap($start,$end,$g_start,$f_start)){
						if (defined $intron_h{$chr}->{$start}){
							my $ori= $intron_h{$chr}->{$start};
							my $post = ";".$pts[0].":".$g_start."-".$f_start;
							$intron_h{$chr}->{$start} = $ori.$post;
						}
						else{
							$intron_h{$chr}->{$start} = $pts[0].":".$g_start."-".$f_start;
						}
					}
				}	
			}
			
			if ($pts[2] eq 'protein'){
				if ($f_start > $m_start){
					print STDERR $lines[$j];
					if( overlap ($start,$end,$m_start,$f_start)){
						if (defined $intron_h{$chr}->{$start}){
							my $ori= $intron_h{$chr}->{$start};
							my $post = ";".$pts[0].":".$m_start."-".$f_start;
							$intron_h{$chr}->{$start} = $ori.$post;
						}
						else{
							$intron_h{$chr}->{$start} = $pts[0].":".$m_start."-".$f_start;
						}
					}
				}
			}
		
			if ($pts[2] eq 'CDS'){
				$last_cor = $f_end;
			}
			
			if ($pts[2] eq 'exon'){
				
				if ($last_cor != 0){
					if (overlap($start,$end,$last_cor,$f_start)){
						
						if ($debug){
							if($i <5){
								print STDERR "\n\nplus:\n",$lines[$j - 1],$lines[$j],"$start,$end,$last_cor,$f_start\n";
							}
						}
						
						if (defined $intron_h{$chr}->{$start}){
							my $ori= $intron_h{$chr}->{$start};
							my $post = ";".$pts[0].":".$last_cor."-".$f_start;
							$intron_h{$chr}->{$start} = $ori.$post;
						}
						else{
							$intron_h{$chr}->{$start} = $pts[0].":".$last_cor."-".$f_start;
						}
					}
					$last_cor = 0;
				}
				
				else{
					if ($g_start == $m_start){
						if ($f_start > $m_start){
							if(overlap($start,$end, $m_start, $f_start)){
								if (defined $intron_h{$chr}->{$start}){
									my $ori= $intron_h{$chr}->{$start};
									my $post = ";".$pts[0].":".$m_start."-".$f_start;
									$intron_h{$chr}->{$start} = $ori.$post;
								}
								else{
									$intron_h{$chr}->{$start} = $pts[0].":".$m_start."-".$f_start;
								}
							}
						}
					}
					else{
						print STDERR "gene_start != mRNA_start\n$g_start != $m_start\n";
					}
				}
			}
			
			if ($j == $#lines){
				if (($f_end < $g_end) and overlap($start,$end,$f_end,$g_end)){
					
					if (defined $intron_h{$chr}->{$start}){
						my $ori= $intron_h{$chr}->{$start};
						my $post = ";".$pts[0].":".$f_end."-".$g_end;
						$intron_h{$chr}->{$start} = $ori.$post;
					}
					else{
						$intron_h{$chr}->{$start} = $pts[0].":".$f_end."-".$g_end;
					}
				}
			}
		}
	}
	
	elsif ($first_pts[6] eq '-'){
		for (my $j = 0; $j <= $#lines;$j++){	
			my $this = $lines[$j];
			next unless ($this =~ /exon|CDS/);
			chomp $this;
			my @pts = split /\t/,$this;
			my ($f_start, $f_end) = @pts[3..4];
			
			
			
			if ($pts[2] eq 'exon'){
				$last_cor = $f_start;
			}
			elsif ($pts[2] eq 'CDS'){
				if ($last_cor != 0){
					if (overlap($start,$end,$f_end,$last_cor)){
						if($debug){
							if($i<10){
								print STDERR "\n\nminus:\n",$lines[$j - 1],$lines[$j],"$start,$end,$f_end,$last_cor\n";
							}
						}
						if (defined $intron_h{$chr}->{$start} ){
							my $ori = $intron_h{$chr}->{$start};
							my $post = ";".$pts[0].":".$f_end."-".$last_cor;
							$intron_h{$chr}->{$start} = $ori.$post;
						}
						else{
							$intron_h{$chr}->{$start} = $pts[0].":".$f_end."-".$last_cor;
						}
					}
					$last_cor = 0;
				}
			}
		}
	}
}

print join("\t",("Chr", "Start", "End", "trim_mean", "Gene", "GeneType", "GeneAnnot",
                  "TE", "TEFamily", "Intergenic","5UTR","exon","intron","3UTR")), "\n";

foreach my $ch (sort keys %peaks_h){
	foreach my $start(sort keys %{$peaks_h{$ch}}){
		my @a = @{$peaks_h{$ch}->{$start}};
		my($five,$exon,$intron,$three) = ("NONE","NONE","NONE","NONE");
		if (defined $five_h{$ch}->{$start} ){
			$five = $five_h{$ch}->{$start} 
		}
		if (defined $cds_h{$ch}->{$start} ){
			$exon = $cds_h{$ch}->{$start} 
		}
		if (defined $intron_h{$ch}->{$start} ){
			$intron = $intron_h{$ch}->{$start} 
		}
		if (defined $three_h{$ch}->{$start} ){
			$three = $three_h{$ch}->{$start} 
		}
		print join ("\t", (@a,$five,$exon,$intron,$three)), "\n";
	}
}

exit;

sub overlap
{
	my ($acor,$bcor,$ccor,$dcor) = @_;
	
	if ($acor>=$ccor && $acor<=$dcor)
	{return 1;}
	if ($bcor>=$ccor && $bcor<=$dcor)
	{return 1;}
	if ($ccor>=$acor && $ccor<=$bcor)
	{return 1;}
	if ($dcor>=$acor && $dcor<=$bcor)
	{return 1;}
	return 0;
}
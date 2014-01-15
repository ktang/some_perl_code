#!/usr/bin/perl -w
=head2
Dec 3,2010(debugging)
Kai TANG
 
 4 NONE the resean is that the gene is TE
version 4
 
 from version 3, version 3 is for list of genes only, this 
 version will be used for output of all the peak region, first use 
 Dr. Liu's add_annotation.pl
 
 
 
 this version is from version 1, not 2
 
 when analysis intron, use only mRNA start and end
 and exon, do not consider CDS.
 
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
my $debug1 = 0;
my $usage = "$0 <input_annotation_file> <gff_file>";

die $usage unless (@ARGV == 2);

my ($inFile,$gffFile) = @ARGV[0..1];

open (IN,$inFile)
	or die "cannot open $inFile:$!";
my %peaks_h;#whole line information of input file;
#my %ends_h;
my %genes_h;# key = gene_name. value = [chr peak_start peak_end]
my %arr_h;# key = gene_name; value is a array contains all the lines assotiated the gene.

my %redundancy_h;#hash store genes which are in at least two peaks.

my %intron_h;
my %cds_h;
my %three_h;
my %five_h;


while (<IN>){
	chomp;
	my @pts = split /\t/;
	$peaks_h{$pts[0]}->{$pts[1]}= \@pts;
	#	$ends_h{$pts[0]}->{$pts[1]} = $pts[2];
	
#	$intron_h{$pts[0]}->{$pts[1]} = 0;
#	$cds_h{$pts[0]}->{$pts[1]} = 0;
#	$three_h{$pts[0]}->{$pts[1]} = 0;
#	$five_h{$pts[0]}->{$pts[1]} = 0;
	if ($pts[4] eq 'NONE'){
		next;
	}	

	if ($pts[4] =~ /;/){
		
		my ($gene1,$gene2) = split /;/,$pts[4];
		
		if(defined $genes_h{$gene1}){
			
			$genes_h{$gene1}->{$redundancy_h{$gene1}} = [$pts[0],$pts[1],$pts[2]];
			++$redundancy_h{$gene1};
		}
			
			
		else{
			$redundancy_h{$gene1} = 1;
				
			$genes_h{$gene1}->{0} = [$pts[0],$pts[1],$pts[2]];
				
		}
		
		if(defined $genes_h{$gene2}){
			$genes_h{$gene2}->{$redundancy_h{$gene2}} = [$pts[0],$pts[1],$pts[2]];
			++$redundancy_h{$gene2};

		}
		else{
			
			$genes_h{$gene2}->{0} = [$pts[0],$pts[1],$pts[2]];
			$redundancy_h{$gene2} = 1;
		}
	}
	
	elsif ($pts[4] ne "NONE" ){
		if(defined $genes_h{$pts[4]}){
			$genes_h{$pts[4]}->{$redundancy_h{$pts[4]}} = [$pts[0],$pts[1],$pts[2]];
			++$redundancy_h{$pts[4]};
		}
		else{
			$genes_h{$pts[4]}->{0} = [$pts[0],$pts[1],$pts[2]];
			$redundancy_h{$pts[4]} = 1;
		}
	}
	
	else{
		die "not a gene nor a NONE,$pts[4]";
	}
}

close(IN);

#my %genes_cor_h;#sotre the start and end of all genes

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
#	$genes_cor_h{$key} = [$pts[3],$pts[4]];
	my $line;
IF:	if (exists $genes_h{$key}){
		
		$line = <GFF>;
		my @mRNAs = split /\t/,$line;
	
		if ( $mRNAs[2] =~ /mRNA|pseudogenic_transcript/){
			push @{$arr_h{$key}},$line;
		}
		else{
			print STDERR "\n\n\nnot mRNA or pseudogenic_transcriptafter gene\n\n$line";
#			die;
		}
	
		while ($line = <GFF>){
			push @{$arr_h{$key}},$line  unless ($line =~ /gene|mRNA/);
			last if ($line =~ /gene|mRNA/);
		}
		
		if ($line =~ /gene/){ 
			my @temp = split /\t/,$line;
			if ($temp[8] =~ /^ID=(AT.G.....);/){
				$key = $1;
				goto IF;
			}
		}

	}		
}
my $i = 0;

foreach my $key (sort keys %genes_h){
	foreach my $num (sort {$a<=>$b}keys %{$genes_h{$key}}){
		$i ++;
		my ($mRNA_start, $mRNA_end) = (0,0);
	
		my @lines = @{$arr_h{$key}};
	
		my ($chr,$start,$end) = @{$genes_h{$key}->{$num}};
#		my ($g_start,$g_end) = @{$genes_cor_h{$key}};
	
		if ($debug){
			if ($i <= 3){
				print STDERR "$i:\n",@lines,"\n";
			}
		}
	
#		my $mRNA_flag = 0;
	
		for (my $j = 0; $j <= $#lines;$j++){
			my $this = $lines[$j];
			chomp $this;
			my @pts = split /\t/,$this;
			
			
			
			if (overlap ($start,$end, $pts[3],$pts[4])){
				
				if ($pts[2] eq "CDS" or $pts[2] eq "pseudogenic_exon"){
					
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
	
#		my ($last_cor, $curr_cor ) = (0,0);
		if($first_pts[2] =~ /mRNA|pseudogenic_transcript/){
			($mRNA_start, $mRNA_end) = ($first_pts[3],$first_pts[4]);
		}
		else{
			print STDERR "\n\n\n\nfirst line is not mRNA pseudogenic_transcript or \n\n\n\n\n$lines[0]";
#			die;
		}

#		my $is_first_exon = 0;
	
		my $last_cor = 0;
	
		if ($first_pts[6] eq '+'){
	
			for (my $j = 1; $j <= $#lines;$j++){
				my $this = $lines[$j];
				next unless ($this =~ /exon/);
				chomp $this;
				my @pts = split /\t/,$this;
				my ($f_start, $f_end) = @pts[3..4];

				if ($pts[2] =~ /exon/){
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
						$last_cor = $f_end;
					}
				
					elsif ($last_cor == 0){ #the first exon
					
						if($f_start > $mRNA_start) {
							print STDERR "\n\nintron is first\n$lines[$j]\n\n";
							if (overlap($start,$end,$mRNA_start,$f_start)){
#								there is an intron before first exon
								$intron_h{$chr}->{$start} = $pts[0].":".$mRNA_start."-".$f_start;
							}
						
						}
						$last_cor = $f_end;
					}
=head2					
					if ($j == $#lines - 1 or $j == $#lines - 2 ){
						if ($f_end < $mRNA_end) {
						
							print STDERR "intron at last\n\n\n$lines[$j]";
						
						
#							if (overlap($start,$end,$f_end,$mRNA_end)){
#								if (defined $intron_h{$chr}->{$start}){
#									my $ori= $intron_h{$chr}->{$start};
#									my $post = ";".$pts[0].":".$f_end."-".$mRNA_end;
#									$intron_h{$chr}->{$start} = $ori.$post;
#								}
#								else{
#									$intron_h{$chr}->{$start} = $pts[0].":".$f_end."-".$mRNA_end;
#								}
#							}
						}
					} 
=cut				
				}
			}
		}
	
		elsif ($first_pts[6] eq '-'){
		
		for (my $j = 1; $j <= $#lines;$j++){	
			my $this = $lines[$j];
			next unless ($this =~ /exon/);
			chomp $this;
			my @pts = split /\t/,$this;
			my ($f_start, $f_end) = @pts[3..4];
			
			if ($pts[2] =~ /exon/){
				
				
			
				if ($last_cor != 0){
					
					if($debug1){
							if( $lines[$j] =~ /AT2G15420/){
								print STDERR "\n\noutside overlap:\n",$lines[$j],"$start,$end,$f_end,$last_cor\n";
							}
						}
					
					if (overlap($start,$end,$f_end,$last_cor)){
						if($debug1){
							if( $lines[$j] =~ /AT2G15420/){
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
					$last_cor = $f_start;
				}
				else{#the first exon
					if ($f_end < $mRNA_end){
						
						print STDERR "minus:intron first:\n\n\n$lines[$j]";
						if (overlap($start, $end, $f_end, $mRNA_end)){
							if (defined $intron_h{$chr}->{$start}){
								my $ori= $intron_h{$chr}->{$start};
								my $post = ";".$pts[0].":".$f_end."-".$mRNA_end;
								$intron_h{$chr}->{$start} = $ori.$post;
							}
							else{
								$intron_h{$chr}->{$start} = $pts[0].":".$f_end."-".$mRNA_end;
							}
						}
					}
						
					$last_cor = $f_start;
				}
				
				if ($j == $#lines){#last exon
					if($mRNA_start < $f_start){
						print STDERR "minus:intron last:\n\n\n$lines[$j]";	
						if ( overlap($start,$end,$mRNA_start, $f_start) ){
							if (defined $intron_h{$chr}->{$start}){
								my $ori= $intron_h{$chr}->{$start};
								my $post = ";".$pts[0].":".$mRNA_start."-".$f_start;
								$intron_h{$chr}->{$start} = $ori.$post;
							}
							else{
								$intron_h{$chr}->{$start} = $pts[0].":".$mRNA_start."-".$f_start;
							}	
						}
					}
				}
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


#if ($debug1){
#	for my $key (sort keys %genes_h){
#		print STDERR "$key\n";
#	}
#}
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

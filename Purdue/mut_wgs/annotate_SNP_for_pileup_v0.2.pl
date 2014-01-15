#!/usr/bin/perl -w
#v0.1
# change how to select snp

#v0.2
# add indel 
=head2
after v2
Feb 22, 2011

extact three letter seq to see 
whether the aa changes.

input format(old)
#Chr	pos	ref	depth	A	C	G	T	SNP_base

new
Chr1    433277  A       6       GGgggg  DIHDDD
Copyright (C) 2010 Kai Tang
version 2 this time add a function to 
determine which position it is in codon(1,2 or 3) 
Nov 23,2010
=cut

use strict;
use Bio::SeqIO;

my $debug = 0;
my $usage = "$0 \n <input_SNP_list> <Col0|C24> <output> STDERR \n\n";
die $usage unless(@ARGV == 3);
my $seq_file = "";

if ($ARGV[1] eq 'Col0'){
	$seq_file = "/Users/tang58/DataBase/TAIR_Col0_genome/5Chr_only_TAIR9_Col0.fas";
}elsif($ARGV[1] eq 'C24'){
	$seq_file = "/Users/tang58/DataBase/C24/C24_TAIR9_5Chr.fas";
}else{
	die $usage;
}

my $snp_file = $ARGV[0];
#my $gff = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_GFF3_CDS.gff";
my $gff = "/Users/tang58/DataBase/TAIR10/GFF/TAIR10_GFF3_CDS.gff";
#my %snps_h;
my %gene_h;
my %seqs;

my $output = $ARGV[2];
die if(-e $output);

open(OUT, ">>$output") or die;

#my %ref_h;#store the best and second best snp base;
my %codon_h;#which position in a codon,1,2 or 3
#my %counts_h;#store the total count{0},best_count{1};second_best{2}
open (SNP,$snp_file)
	or die "cannot open $snp_file:$!";
open (CDS,$gff)
	or die "cannot open $gff:$!";

my $seq_in = Bio::SeqIO->new('-file' => "$seq_file",
                             '-format' => 'fasta');

while (my $inseq = $seq_in->next_seq){
	my $chr = lc $inseq->id;
	$seqs{$chr} = $inseq->seq;
}

#my $seq_debug = substr( $seqs{chr1}, 7156, 10);
#print STDERR $seq_debug,"\n"; 


my @SNP = <CDS>;
close(CDS);
my %snps;
my %records;
my %indels;

while (<SNP>){
	chomp;
	my @pts = split /\t/;
	my $chr = lc $pts[0];
	my $pos = $pts[1];
	my $dep = $pts[3];
	my $seq = $pts[4];
	#$snps{$chr}->{$pts[1]} = [@pts[0..2], $pts[8]];
	my $snp_base = "N";
	
	my $star_num = ( $seq =~ tr/\*/\*/ );
	
	if($star_num / $dep >= 0.3){
		$indels{$chr}->{$pos} = 1;
		$records{$chr}->{$pos}  = $_;
	}
	
#	if (  ((( $seq =~ tr/a/a/) + ($seq =~ tr/A/A/) ) / $dep ) >= 0.5 ){
#		$snp_base = "A";
#	}elsif(((( $seq =~ tr/c/c/) + ($seq =~ tr/C/C/) ) / $dep ) >= 0.5 ){
#		$snp_base = "C";
#	}elsif(((( $seq =~ tr/g/g/) + ($seq =~ tr/G/G/) ) / $dep ) >= 0.5 ){
#		$snp_base = "G";
#	}elsif(((( $seq =~ tr/t/t/) + ($seq =~ tr/T/T/) ) / $dep ) >= 0.5 ){
#		$snp_base = "T";
#	}else{
#		print  $_, "\n";
#	}

	#v0.1
	my $A_num = ( $seq =~ tr/Aa/Aa/  );
	my $C_num = ( $seq =~ tr/Cc/Cc/  );
	my $G_num = ( $seq =~ tr/Gg/Gg/  );
	my $T_num = ( $seq =~ tr/Tt/Tt/  );
	
	my $max = max($A_num, $C_num, $G_num, $T_num);
	my $hit_num = 0;
	if( $max == $A_num) {$hit_num ++; $snp_base = "A"; }
	if( $max == $C_num) {$hit_num ++; $snp_base = "C"; }
	if( $max == $G_num) {$hit_num ++; $snp_base = "G"; }
	if( $max == $T_num) {$hit_num ++; $snp_base = "T"; }
	
	unless ($hit_num == 1){
		print STDERR join("\t", ($max, $A_num, $C_num, $G_num, $T_num, $_)), "\n" if($max >= 3);
		$snp_base = "N";
	}
	
	if($snp_base ne "N" and $max >= 3){
		$snps{$chr}->{$pos} = [@pts[0..2], $snp_base];
		$records{$chr}->{$pos}  = $_;
	}
#	if ($pts[2] ne $pts[5] or $pts[12] >= 5){
#		$ref_h{$chr}->{$pts[1]} = $pts[2];
#		$snps_h{$chr}->{$pts[1]}->{1}= $pts[5];
#		$snps_h{$chr}->{$pts[1]}->{2}= $pts[9];
#		$counts_h{$chr}->{$pts[1]}-> {0} = $pts[13];
#		$counts_h{$chr}->{$pts[1]}-> {1} = $pts[8];
#		$counts_h{$chr}->{$pts[1]}-> {2} = $pts[12];
#	}
}
close(SNP);

for (my $i = 0; $i<=$#SNP; $i++){
	
	my $this = $SNP[$i];
	chomp $this;
	my @temp = split /\t/, $this;
	my $chr = lc $temp[0];
	my ($start, $end) = ( $temp[3],  $temp[4]);
	
	foreach my $pos (keys %{$snps{$chr}}){
		if ($pos >= $temp[3] and $pos <= $temp[4]){
				my $snp = ${$snps{$chr}->{$pos}}[3];
				
				
				my $codon_seq = "";
				my $codon_sub_seq = "";
	     		my @parents = split /,/,$temp[8];
				my @genes = split /=/, $parents[0];
				$gene_h{$chr}->{$pos} = $genes[1];
				if ($temp[6] eq "+"){
					my $codon_pos = $codon_h{$chr}->{$pos} = ($pos - $temp[3] + 3 - $temp[7]) % 3 + 1;
					
					
					
					#extract codon seq
					if($codon_pos == 1){
						
						if($end - $pos >= 2){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 1, 3);
						}elsif($end - $pos == 1){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 2);
							my @a = split "\t", $SNP[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base;
						}elsif($end - $pos == 0){#next two base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $SNP[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 2);
							$codon_seq .= $next_base;
						}else{
							die "codon == 1\n";
						}			
						
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 0, 1, $snp);
								
					}elsif($codon_pos == 2){
						
						if($pos > $start and $pos < $end){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 2, 3);
						}elsif($pos == $start and $pos < $end){#first base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 2);
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $end and $pos > $start){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							my @a = split "\t", $SNP[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base;
						}else{
							die "codon == 2\n";
						}
						
					#	$codon_sub_seq = substr($codon_seq, 1, 1, $snp);
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 1, 1, $snp);
						
					}elsif($codon_pos == 3){
						
						if($pos - $start >= 2){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 3, 3);
						}elsif($pos - $start == 1){# first base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							#my $next_cds = $SNP[$i + 1];
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $start){# first two base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[4] - 2, 2);
							$codon_seq = $next_base.$codon_seq ;
						}else{
							die "codon == 3\n";
						}
						
						#$codon_sub_seq = substr($codon_seq, 2, 1, $snp);	
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 2, 1, $snp);
						
					}else{
						die "wrong codon_position\n";
					}
				
					my $codon_ori = Bio::Seq->new(-seq => $codon_seq);
					my $codon_sub = Bio::Seq->new(-seq => $codon_sub_seq);
#					$prot_obj = $dna_obj->translate
					my $pro_ori = $codon_ori->translate;
					my $aa_ori = $pro_ori->seq;
					my $pro_sub = $codon_sub->translate;
					my $aa_sub = $pro_sub->seq;
					if ($aa_ori ne $aa_sub){
						print OUT join("\t",($chr,$pos,$genes[1],"+",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub, $records{$chr}->{$pos})), "\n";
					}
#					if($debug){
#							print STDERR join("\t",($chr,$pos,"+",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub,$snp)), "\n";
#					}

				}
				
				elsif ($temp[6] eq "-"){
					my $codon_pos = $codon_h{$chr}->{$pos} = ($temp[4] - $pos + 3 - $temp[7]) % 3 + 1;
					
#					print STDERR "$chr,$pos,$codon_pos\n";
					
					#extract codon seq
					if($codon_pos == 1){
						
						if($pos - $start >= 2){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 3, 3);
						}elsif($pos - $start == 1){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							
							my @a = split "\t", $SNP[$i + 1];# next CDS
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							# notic easy to make a mistake as $codon_seq .= $next_base
							$codon_seq = $next_base.$codon_seq;
						}elsif($pos == $start){#next two base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $SNP[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[4] - 2, 2);
							$codon_seq = $next_base.$codon_seq;
						}else{
							die "codon == 1\n";
						}			
						
#						$codon_sub_seq = substr($codon_seq, 2, 1, $snp);
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 2, 1, $snp);								
					}elsif($codon_pos == 2){
						
						if($pos > $start and $pos < $end){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 2, 3);
						}elsif($pos == $start and $pos < $end){#first base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 2);
							my @a = split "\t", $SNP[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $end and $pos > $start){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base;
						}else{
							die "codon == 2\n";
						}
						
					#	$codon_sub_seq = substr($codon_seq, 1, 1, $snp);
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 1, 1, $snp);
						
					}elsif($codon_pos == 3){#321 in plus strand
						
						if($end - $pos >= 2){#normal
							$codon_seq = substr($seqs{$chr}, $pos - 1, 3);
						}elsif($end - $pos == 1){# first base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 2);
							#my $next_cds = $SNP[$i + 1];
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base ;
						}elsif($pos == $end){# first two base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $SNP[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[3] - 1, 2);
							$codon_seq .= $next_base ;
						}else{
							die "codon == 3\n";
						}
						
					#	$codon_sub_seq = substr($codon_seq, 0, 1, $snp);
						$codon_sub_seq = $codon_seq;
						substr($codon_sub_seq, 0, 1, $snp);
						
					}else{
						die "wrong codon_position\n";
					}

					my $codon_ori = Bio::Seq->new(-seq => $codon_seq);
					my $codon_sub = Bio::Seq->new(-seq => $codon_sub_seq);
					my $rev_ori = $codon_ori->revcom;
					my $rev_sub = $codon_sub->revcom;
#					$prot_obj = $dna_obj->translate
					my $pro_ori = $rev_ori->translate;
					my $aa_ori = $pro_ori->seq;
					my $pro_sub = $rev_sub->translate;
					my $aa_sub = $pro_sub->seq;
					if ($aa_ori ne $aa_sub){
						print OUT join("\t",($chr,$pos,$genes[1],"-",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub, $records{$chr}->{$pos})), "\n";
					}
					
#					if($debug){
#							print STDERR join("\t",($chr,$pos,"-",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub,$snp)), "\n";
#					}
					
				}
				else{die "neigth + nor -"}
			
		}
	}
	
}

for (my $i = 0; $i<=$#SNP; $i++){
	
	my $this = $SNP[$i];
	chomp $this;
	my @temp = split /\t/, $this;
	my $chr = lc $temp[0];
	my ($start, $end) = ( $temp[3],  $temp[4]);
	foreach my $pos (keys %{$indels{$chr}}){
		if ($pos >= $temp[3] and $pos <= $temp[4]){
			my @parents = split /,/,$temp[8];
			my @genes = split /=/, $parents[0];
			$gene_h{$chr}->{$pos} = $genes[1];
			
			print OUT join("\t",($chr,$pos,$genes[1], $records{$chr}->{$pos})), "\n";
		
		}
	}
}

exit;

sub max{
	my @a = @_;
	my @b = sort {$b <=> $a} @a;
	return $b[0];
}


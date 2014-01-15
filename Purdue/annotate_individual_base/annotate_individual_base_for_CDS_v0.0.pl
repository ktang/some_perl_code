#!/usr/bin/perl -w

#v0.0
# specific for Xiaohong's project
# input format

#chr     pos     ref     mut_WT_0        mut_WT_150      mut_ddc_0       mut_ddc_150     mut_nrpe1_0     mut_nrpe1_150   seq_WT_0        seq_WT_150      seq_ddc_0       seq_ddc_150     seq_nrpe1_0     seq_nrpe1_150   type    type_code
#0	 1	 2	 3		 4			5		6		7		8		9		10		11		12		13		14	 15		16


###########################
# modified from /Users/tang58/scripts_all/perl_code/Purdue/mut_wgs/annotate_SNP_for_pileup_v0.2.pl
############################
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

#input format

my $debug = 0;
my $usage = "$0 \n <input_SNP_list> <Col0|C24> STDOUT \n\n";
die $usage unless(@ARGV == 2);
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

#my $output = $ARGV[2];
#die if(-e $output);
#open(OUT, ">>$output") or die;

#my %ref_h;#store the best and second best snp base;
my %codon_h;#which position in a codon,1,2 or 3
#my %counts_h;#store the total count{0},best_count{1};second_best{2}
open (SNP, $snp_file)	or die "cannot open $snp_file:$!";
open (CDS, $gff)	or die "cannot open $gff:$!";


my $seq_in = Bio::SeqIO->new('-file' => "$seq_file",
                             '-format' => 'fasta');

while (my $inseq = $seq_in->next_seq){
	my $chr = lc $inseq->id;
	$chr =~ s/chr//i;
	$seqs{$chr} = $inseq->seq;
}


my @CDS = <CDS>;
close(CDS);
my %snps;
my %records;
my %indels;

my $head = <SNP>;

while (<SNP>){
	chomp;
	my @pts = split /\t/;
	my $chr = lc $pts[0];
	$chr =~ s/chr//i;
	my $pos = $pts[1];
#	my $dep = $pts[3];
#	my $seq = $pts[4];
	#$snps{$chr}->{$pts[1]} = [@pts[0..2], $pts[8]];
	my $snp_base = "N";
	
	for my $i (3..8){
		next if ( $pts[$i] eq "NONE");
		if( $pts[$i] =~ /\S=>([ACTG])/ ){
			my $tmp = $1;
			print STDERR "mult:", $_, "\n" if ( $snp_base ne "N" and  $snp_base ne $tmp );
			$snp_base = $tmp;
		}
	}
	
	if( $snp_base eq "N" ){
		$snp_base = $pts[2];
	}
	
	$snps{$chr}->{$pos} = [@pts[0..2], $snp_base];
	$records{$chr}->{$pos}  = $_;
	

}
close(SNP);

# @CDS is list for CDS
for (my $i = 0; $i<=$#CDS; $i++){
	
	my $this = $CDS[$i];
	chomp $this;
	my @temp = split /\t/, $this;
	my $chr = lc $temp[0];
	$chr =~ s/chr//i;
	my ($start, $end) = ( $temp[3],  $temp[4]);

	foreach my $pos (keys %{$snps{$chr}}){
		if ($pos >= $temp[3] and $pos <= $temp[4]){ # SNP is in one of the CDS
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
							my @a = split "\t", $CDS[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base;
						}elsif($end - $pos == 0){#next two base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $CDS[$i + 1];
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
							my @a = split "\t", $CDS[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $end and $pos > $start){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							my @a = split "\t", $CDS[$i + 1];
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
							#my $next_cds = $CDS[$i + 1];
							my @a = split "\t", $CDS[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $start){# first two base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $CDS[$i - 1];
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
					my $pro_ori = $codon_ori->translate;
					my $aa_ori = $pro_ori->seq;
					my $pro_sub = $codon_sub->translate;
					my $aa_sub = $pro_sub->seq;
					if ($aa_ori ne $aa_sub){
						print join("\t",($chr,$pos,$genes[1],"+",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub, $records{$chr}->{$pos})), "\n";
					}
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
							
							my @a = split "\t", $CDS[$i + 1];# next CDS
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							# notic easy to make a mistake as $codon_seq .= $next_base
							$codon_seq = $next_base.$codon_seq;
						}elsif($pos == $start){#next two base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $CDS[$i + 1];
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
							my @a = split "\t", $CDS[$i + 1];
							my $next_base = substr($seqs{$chr}, $a[4] -1, 1);
							$codon_seq = $next_base.$codon_seq ;
						}elsif($pos == $end and $pos > $start){#last base in next cds
							$codon_seq = substr($seqs{$chr}, $pos - 2, 2);
							my @a = split "\t", $CDS[$i - 1];
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
							#my $next_cds = $CDS[$i + 1];
							my @a = split "\t", $CDS[$i - 1];
							my $next_base = substr($seqs{$chr}, $a[3] -1, 1);
							$codon_seq .= $next_base ;
						}elsif($pos == $end){# first two base in former cds
							$codon_seq = substr($seqs{$chr}, $pos - 1, 1);
							my @a = split "\t", $CDS[$i - 1];
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
					my $pro_ori = $rev_ori->translate;
					my $aa_ori = $pro_ori->seq;
					my $pro_sub = $rev_sub->translate;
					my $aa_sub = $pro_sub->seq;
					if ($aa_ori ne $aa_sub){
						print join("\t",($chr,$pos,$genes[1],"-",$codon_seq,$aa_ori,$codon_sub_seq,$aa_sub, $records{$chr}->{$pos})), "\n";
					}
				}
				else{
					die "neigth + nor -"
				}
		}
	}
}

# no indel this time
=head
for (my $i = 0; $i<=$#CDS; $i++){
	
	my $this = $CDS[$i];
	chomp $this;
	my @temp = split /\t/, $this;
	my $chr = lc $temp[0];
	my ($start, $end) = ( $temp[3],  $temp[4]);
	foreach my $pos (keys %{$indels{$chr}}){
		if ($pos >= $temp[3] and $pos <= $temp[4]){
			my @parents = split /,/,$temp[8];
			my @genes = split /=/, $parents[0];
			$gene_h{$chr}->{$pos} = $genes[1];
			
			print STDERR join("\t",($chr,$pos,$genes[1], $records{$chr}->{$pos})), "\n";
		
		}
	}
}
=cut

exit;

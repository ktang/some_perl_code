#!/usr/bin/perl -w

=head2

all C
percentage means
(rdd_mC_plus + rdd_mC_minus - Col_mC_plus - Col_mC_minus)/(TAIR9_C_plus + TAIR9_C_minus)



v2.0_score
not number of mCG, but use score
v2.0
only number, no strand specific,


this is verson 1.9
strand for genes.

v1.8
focus on CG only

output format
no score output
average of forward and reverse strand

#Chr |Start |end |+_strand_C_number |rdd_mC_No_of_+ |rdd_mC_score_of_+ |col0_mC_No_of_+ |col0_mC_score_of_+ |-_strand_C_number |rdd_mC_No_of_- |rdd_mC_score_of_- |col0_mC_No_of_- col0_mC_score_of_-
number is not int, can be persentage

use bed file as input

output format should be
chr	.	.	start	end	total_C	mC_rdd	mc_col	diff_percentage

the difference from verson 4 is that 
this script counts the average of minus strand and plus strand

in Ecker data.

exclude region which is not sequenced.

this script is use the TAIR9_chr.fasta file to
calculate the CG,CHG and CHH number within peak region,
then use Ecker data to calculate methylated C in two lines
(WT and rdd); then use a cutoff to select the difference of 
the two. (rdd - WT >= cutoff) 

/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_col0_tair9_no_warn.gff

notice: the bed file must be sorted
=cut

use Bio::SeqIO;
use strict;

my $usage = "$0 <peak_bed_file> ";

die $usage unless (@ARGV == 1 ); 

my $bed_file = $ARGV[0];



my $chr_file = "/Users/tang58/DataBase/TAIR_Col0_genome/5Chr_only_TAIR9_Col0.fas";

my $ecker_col0 = "/Users/tang58/Nimblegen_project/Ecker_data/mc_col0_tair9_no_warn.gff";

my $ecker_rdd = "/Users/tang58/Nimblegen_project/Ecker_data/mc_rdd_tair9_no_warn.gff";

my $chr_in = Bio::SeqIO -> new ('-file' => $chr_file,
							    '-format' => "fasta");

open (BED, "<$bed_file")
	or die "cannot open $bed_file :$!";
	
open(WT, "<$ecker_col0")
	or die "cannot open $ecker_col0: $!";
	
open (RDD, "<$ecker_rdd")
	or die "cannot open $ecker_rdd: $!";
	
my @beds = <BED>;
close (BED);

my @col0s = <WT>;
close(WT);

my @rdds = <RDD>;
close(RDD);

my %chrs;
while (my $chr_seq_obj = $chr_in -> next_seq)
{
	my $chr_name = lc($chr_seq_obj -> id);
	my $chr_seq = 	$chr_seq_obj -> seq;
	$chrs{$chr_name} = $chr_seq;
}

$chr_in->close;

my %counts;
my %counts_minus;

for (my $i = 0; $i <= $#beds; $i++){
	my $thisPeak = $beds[$i];
	chomp $thisPeak;
	my @pts_peak = split "\t", $thisPeak;
	my $peak_chr = lc($pts_peak[0]);
	my $peak_start = $pts_peak[1] + 0;
	my $peak_end = $pts_peak[2] + 0;
	my $peak_length = $peak_end - $peak_start + 1;
	
	my $peak_seq = substr($chrs{$peak_chr},$peak_start - 1, $peak_length);
	my $plus = ($peak_seq =~ tr/C/C/);
	my $minus = ($peak_seq =~ tr/G/G/);
	if ($plus == 0) {$plus = 0}
	if ($minus == 0) {$minus = 0}	
	$counts{$peak_chr}->{$pts_peak[1]} = $plus;
	$counts_minus{$peak_chr}->{$pts_peak[1]} = $minus;
}


###########################################################################
my ($b,$c) = (0,0);

my %col_met;
my %col_met_score;
my %col_met_minus;
my %col_met_minus_score;

LOOP: while ( $b <= $#beds and $c <= $#col0s )
{
	my $thisPeak = $beds[$b];
	my $thisC = $col0s[$c];
	chomp $thisPeak;
	chomp $thisC;
	my @pts_bed = split "\t",$thisPeak;
	my @pts_C = split "\t", $thisC;
	my $pos_C = $pts_C[3];
	my $start = $pts_bed[1] + 0;
	my $end = $pts_bed[2] + 0;
	
	my $chr = lc($pts_bed[0]);
	#my $strand = $pts_bed[6];
	if ( lc($pts_bed[0]) eq lc($pts_C[0]))
	{
		if ( $pos_C < $start){
			$c++;
			next LOOP;	
		}
		elsif ( $pos_C > $end){
			$b++;
			next LOOP;	
		}
		elsif ( ( $pos_C >= $start) and ($pos_C <= $end))
		{
						
			$col_met{$chr}->{$pts_bed[1]}  += ($pts_C[6] / 100);
			$col_met{$chr}->{$pts_bed[1]}  += ($pts_C[7] / 100);
			$c++;
			next LOOP;
		}
	}		
	elsif (lc($pts_bed[0]) lt lc($pts_C[0])){
		$b++;
		next LOOP;	
	}
	elsif (lc($pts_bed[0]) gt lc($pts_C[0]) ) {
		$c++;
		next LOOP;	
	}
}

###########################################################################
($b,$c) = (0,0);

my %rdd_met;
my %rdd_met_score;
my %rdd_met_minus;
my %rdd_met_minus_score;
LOOP: while ( $b <= $#beds and $c <= $#rdds )
{
	my $thisPeak = $beds[$b];
	my $thisC = $rdds[$c];
	chomp $thisPeak;
	chomp $thisC;
	my @pts_bed = split "\t",$thisPeak;
	my @pts_C = split "\t", $thisC;
	my $pos_C = $pts_C[3];
	my $start = $pts_bed[1] + 0;
	my $end = $pts_bed[2] + 0;
	
	my $chr = lc($pts_bed[0]);
	if ( lc($pts_bed[0]) eq lc($pts_C[0]))
	{
		if ( $pos_C < $start){
			$c++;
			next LOOP;	
		}
		elsif ( $pos_C > $end){
			$b++;
			next LOOP;	
		}
		elsif ( ( $pos_C >= $start) and ($pos_C <= $end))
		{
			$rdd_met{$chr}->{$pts_bed[1]} += ($pts_C[6] / 100);
			$rdd_met{$chr}->{$pts_bed[1]} += ($pts_C[7] / 100);
			$c++;
			next LOOP;
		}
	}		
	elsif (lc($pts_bed[0]) lt lc($pts_C[0])){
		$b++;
		next LOOP;	
	}
	elsif (lc($pts_bed[0]) gt lc($pts_C[0]) ) {
		$c++;
		next LOOP;	
	}
}

###########################################################################


for (my $i = 0; $i <= $#beds; $i++){
	my $this = $beds[$i];
	chomp $this;
	my @pts =  split "\t", $this;
	my $chr = lc($pts[0]);
	
	
	
	my $C = $counts{$chr}->{$pts[1]} + $counts_minus{$chr}->{$pts[1]};
	
	my ($mC_col0, $mC_rdd)= (0, 0);
	
	if (defined $col_met{$chr}->{$pts[1]}){
		$mC_col0 = $col_met{$chr}->{$pts[1]};
	}
	
	if (defined $rdd_met{$chr}->{$pts[1]}){
		$mC_rdd = $rdd_met{$chr}->{$pts[1]};
	}
	

	my $diff = 0;
	if ($C != 0){
		$diff = ($mC_rdd - $mC_col0) / $C;
	}
	else{
		print STDERR "C = 0\n ", join("\t", ($pts[0], $pts[1], $pts[2], $C, $mC_rdd, $mC_col0, $diff)), "\n";
	}
	print join("\t", ($pts[0], $pts[1], $pts[2], $C, $mC_rdd, $mC_col0, $diff)), "\n";
}
exit;
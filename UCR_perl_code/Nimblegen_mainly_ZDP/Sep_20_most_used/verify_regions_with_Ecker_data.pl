#!/usr/bin/perl -w

=head2
this script is use the TAIR9_chr.fasta file to
calculate the CG,CHG and CHH number within peak region,
then use Ecker data to calculate methylated C in two lines
(WT and rdd); then use a cutoff to select the difference of 
the two. (rdd - WT >= cutoff) 

/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_col0_tair9_no_warn.gff

notice: the bed file must be sorted
=cut

BEGIN{unshift @INC, "/Library/Perl/5.10.0/BioPerl";}

use Bio::SeqIO;
use strict;

my $usage = "$0 <peak_bed_file> <cutoff_cg> <cutoff_chg> <cutoff_chh> <and_or>";

die $usage unless ((@ARGV == 5 ) and (($ARGV[4] eq "and" ) or($ARGV[4] eq "or")));

my $bed_file = $ARGV[0];

my $cutoff_cg = $ARGV[1];
my $cutoff_chg = $ARGV[2];
my $cutoff_chh = $ARGV[3];


my $chr_file = "/Users/kaitang/Desktop/TAIR/chromosomes/TAIR9_chr_all.fas";

my $ccker_col0 = "/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_col0_tair9_no_warn.gff";

my $ccker_rdd = "/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_rdd_tair9_no_warn.gff";

my $chr_in = Bio::SeqIO -> new ('-file' => $chr_file,
			    '-format' => "fasta");

open (BED, "<$bed_file")
	or die "cannot open $bed_file :$!";
	
open(WT, "<$ccker_col0")
	or die "cannot open $ccker_col0: $!";
	
open (RDD, "<$ccker_rdd")
	or die "cannot open $ccker_rdd: $!";
	
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

for (my $i = 0; $i <= $#beds; $i++)
{
	my $thisPeak = $beds[$i];
	chomp $thisPeak;
	my @pts_peak = split "\t", $thisPeak;
	my $peak_chr = lc($pts_peak[0]);
	my $peak_start = $pts_peak[1] + 0;
	my $peak_end = $pts_peak[2] + 0;
	my $peak_length = $peak_end - $peak_start + 1;
	
	my $peak_seq = substr($chrs{$peak_chr},$peak_start - 1, $peak_length+2);
	my @peak_seq_array = split "",$peak_seq;
	
	for (my $j = 0; $j < $peak_length; $j++)
	{
			if ($peak_seq_array[$j] eq "C")
			{
				if ($peak_seq_array[$j+1] eq "G" )
				{$counts{$peak_chr}->{$pts_peak[1]}->{"CG"}++}
				else
				{
					if($peak_seq_array[$j+2] eq "G")
					{
						$counts{$peak_chr}->{$pts_peak[1]}->{"CHG"}++;
					}
					else{$counts{$peak_chr}->{$pts_peak[1]}->{"CHH"}++}
				}
			}
	}
}


###########################################################################
my ($b,$c) = (0,0);

my %col_met;

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
	my $cnd = $pts_bed[2] + 0;
	
	my $chr = lc($pts_bed[0]);
	#my $strand = $pts_bed[6];
	if ( lc($pts_bed[0]) eq lc($pts_C[0]))
	{
		if ( $pos_C < $start){
			$c++;
			next LOOP;	
		}
		elsif ( $pos_C > $cnd){
			$b++;
			next LOOP;	
		}
		elsif ( ( $pos_C >= $start) and ($pos_C <= $cnd))
		{
			if ( $pts_C[6] > 0 )
			{
						if ($pts_C[5] eq "CG") { $col_met{$chr}->{$pts_bed[1]}->{"CG"}++ }
						if ($pts_C[5] eq "CHG") { $col_met{$chr}->{$pts_bed[1]}->{"CHG"}++ }
						if ($pts_C[5] eq "CHH") { $col_met{$chr}->{$pts_bed[1]}->{"CHH"}++ }
			}
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
	my $cnd = $pts_bed[2] + 0;
	
	my $chr = lc($pts_bed[0]);
	#my $strand = $pts_bed[6];
	if ( lc($pts_bed[0]) eq lc($pts_C[0]))
	{
		if ( $pos_C < $start){
			$c++;
			next LOOP;	
		}
		elsif ( $pos_C > $cnd){
			$b++;
			next LOOP;	
		}
		elsif ( ( $pos_C >= $start) and ($pos_C <= $cnd))
		{
			if ( $pts_C[6] > 0 )
			{
						if ($pts_C[5] eq "CG") { $rdd_met{$chr}->{$pts_bed[1]}->{"CG"}++ }
						if ($pts_C[5] eq "CHG") { $rdd_met{$chr}->{$pts_bed[1]}->{"CHG"}++ }
						if ($pts_C[5] eq "CHH") { $rdd_met{$chr}->{$pts_bed[1]}->{"CHH"}++ }
			}
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
for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my @pts =  split "\t", $this;
	my $chr = lc($pts[0]);
	my $cg = $counts{$chr}->{$pts[1]}->{"CG"};
	my $chg = $counts{$chr}->{$pts[1]}->{"CHG"}; 
	my $chh = $counts{$chr}->{$pts[1]}->{"CHH"};
	my ($cg_col0 ,$chg_col0 ,$chh_col0  ) = (0,0,0);
	if ( (exists $col_met{$chr}) and (exists $col_met{$chr}->{$pts[1]}) )
	{
		if (exists $col_met{$chr}->{$pts[1]}->{"CG"} )
		{
			$cg_col0 = $col_met{$chr}->{$pts[1]}->{"CG"};
		}
		
		if (exists $col_met{$chr}->{$pts[1]}->{"CHG"} )	
		{
			$chg_col0 = $col_met{$chr}->{$pts[1]}->{"CHG"};
		}
		if (exists $col_met{$chr}->{$pts[1]}->{"CHH"} )
		{
			$chh_col0 = $col_met{$chr}->{$pts[1]}->{"CHH"};
		}
	}
	
	
	my ($cg_rdd ,$chg_rdd  ,$chh_rdd   ) = (0,0,0);
	if ( (exists $rdd_met{$chr}) and (exists $rdd_met{$chr}->{$pts[1]}) )
	{
		if (exists $rdd_met{$chr}->{$pts[1]}->{"CG"} )
		{
			$cg_rdd = $rdd_met{$chr}->{$pts[1]}->{"CG"};
		}
		
		if (exists $rdd_met{$chr}->{$pts[1]}->{"CHG"} )	
		{
			$chg_rdd = $rdd_met{$chr}->{$pts[1]}->{"CHG"};
		}
		if (exists $rdd_met{$chr}->{$pts[1]}->{"CHH"} )
		{
			$chh_rdd = $rdd_met{$chr}->{$pts[1]}->{"CHH"};
		}
	}
#	if ( $cg_col0 == undef) {$cg_col0 = 0}
#	if ( $chg_col0 == undef) {$chg_col0 = 0}
#	if ( $chh_col0 == undef) {$chh_col0 = 0}
	
#	my $cg_rdd = $rdd_met{$chr}->{$pts[1]}->{"CG"};
#	my $chg_rdd = $rdd_met{$chr}->{$pts[1]}->{"CHG"};
#	my $chh_rdd = $rdd_met{$chr}->{$pts[1]}->{"CHH"};
#	if( $cg_rdd == undef) {$cg_rdd = 0}
#	if( $chg_rdd == undef) {$chg_rdd = 0}
#	if( $chh_rdd == undef) {$chh_rdd = 0}
 
	my $diff_cg = ($cg_rdd - $cg_col0)/$cg;
	my $diff_chg = ($chg_rdd - $chg_col0)/$chg;
	my $diff_chh = ($chh_rdd - $chh_col0)/$chh;
	
#	print STDERR "$cg\t$chg\t$chh\t$cg_col0\t$chg_col0\t$chh_col0\t$cg_rdd\t$chg_rdd\t$chh_rdd\n";
#	print STDERR "$diff_cg\t$diff_chg\t$diff_chh\n";
	
	if ( $ARGV[4] eq "and")
	{
		if ( ($diff_cg >= $cutoff_cg) and ($diff_chg >= $cutoff_chg) and ($diff_chh >= $cutoff_chh))
		{
			print $this,"\n";	
		}
	}
	elsif ($ARGV[4] eq "or")
	{
		if ( ($diff_cg >= $cutoff_cg) or ($diff_chg >= $cutoff_chg) or ($diff_chh >= $cutoff_chh))
		{
			print $this,"\n";	
		}	
	}
}

exit;
#!/usr/bin/perl -w

=head2
this script only count the mC number in the Ecker data.
#the problem is that G in the + strand is C in the 
# - strand
# the fact is that most(1869913/1869944) is mapped to the + strands
# so I should count the C in the - strands, which is G in + strands
# this time use a cutoff to determine whether a C is methylated.

this script is change from verify_regions_with_Ecker_data_v6.pl


use bed(with trimmean in the last column) file as input

output format should be
chr	start	end	trimmean	length	mC	



exclude region which is not sequenced.

this script is use the TAIR9_chr.fasta file to
calculate the CG,CHG and CHH number within peak region,
then use Ecker data to calculate methylated C in two lines
(WT and rdd); then use a cutoff to select the difference of 
the two. (rdd - WT > cutoff) 

/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_col0_tair9_no_warn.gff

notice: the bed file must be sorted
=cut


use strict;

my $usage = "$0 <peak_bed_file> <Ecker_data> <cutoff_0_100>";

die $usage unless (@ARGV == 3 ); #and (($ARGV[4] eq "and" ) or($ARGV[4] eq "or")));

my $bed_file = $ARGV[0];
my $cutoff = $ARGV[2];



#my $chr_file = "/Users/kaitang/Desktop/TAIR/chromosomes/TAIR9_chr_all.fas";

#my $ccker_col0 = "/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_col0_tair9_no_warn.gff";

#my $ccker_rdd = "/Users/kaitang/Desktop/Nimblegen/Ecker_data/mC_data_modified_TAIR8_to_TAIR9/mc_rdd_tair9_no_warn.gff";

my $ecker_file = $ARGV[1];

#my $chr_in = Bio::SeqIO -> new ('-file' => $chr_file,
	#		    '-format' => "fasta");

open (BED, "<$bed_file")
	or die "cannot open $bed_file :$!";
	
open (ECKER, "$ecker_file")
	or die "cannot open $ecker_file:$!";
	
my @beds = <BED>;
close (BED);

my @eckers = <ECKER>;
close (ECKER);

=head2
#my %chrs;
#while (my $chr_seq_obj = $chr_in -> next_seq)
#{
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
	
	my $peak_seq = substr($chrs{$peak_chr},$peak_start - 1, $peak_length);
	#	my @peak_seq_array = split "",$peak_seq;
	#############################################
	#
	#tr/G/G/ not C
	############################################
	$counts{$peak_chr}->{$pts_peak[1]} = ($peak_seq =~ tr/G/G/);
=cut
=head1
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
=cut


###########################################################################
my ($b,$e) = (0,0);

my %mC;

LOOP: while ( $b <= $#beds and $e <= $#eckers )
{
	my $thisPeak = $beds[$b];
	my $thisC = $eckers[$e];
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
			$e++;
			next LOOP;	
		}
		elsif ( $pos_C > $end){
			$b++;
			next LOOP;	
		}
		
		
		############################################################
		elsif ( ( $pos_C >= $start) and ($pos_C <= $end))
		{
			if ( $pts_C[7] > $cutoff )
			{
								$mC{$chr}->{$pts_bed[1]}++;
			}
			$e++;
			next LOOP;
		}
	}		
	elsif (lc($pts_bed[0]) lt lc($pts_C[0])){
		$b++;
		next LOOP;	
	}
	elsif (lc($pts_bed[0]) gt lc($pts_C[0]) ) {
		$e++;
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

#	my $c = $counts{$chr}->{$pts[1]};
	
	my $met_C = 0;
	if ((exists $mC{$chr}) and (exists $mC{$chr}->{$pts[1]}))
	{
		$met_C = $mC{$chr}->{$pts[1]};
	}
	
=head2		
	my $percentage = 0;
	if ($c == 0 )
	{
		if ($met_C != 0)
		{
				print STDERR "$chr\t$pts[1]\t$pts[2]\t$pts[3]n";
				die;
		}
	}
	
	else
	{
	 	 $percentage=$met_C/$c;
	}
=cut
	my $len = $pts[2] - $pts[1] + 1 ;
	
	print "$chr\t$pts[1]\t$pts[2]\t$pts[3]\t$len\t$met_C\n";
	#if ( ($i < 10) or ($i > ($#beds -5)))
	#{
	#	print STDERR "$pts[0]\t$pts[1]\tC=$c\tcol0=$mC_col0\trdd=$mC_rdd\n";	
	#}
}

exit;

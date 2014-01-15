#!/usr/bin/perl -w

use strict;

#START HERE MAIN PROGRAM

my $usage = "$0 <Ecker_data> <peaks_gff>";
die $usage unless (@ARGV == 2);

my $ecker_file = $ARGV[0];
my $peakfilenm=$ARGV[1];
my %nums;

#load gff lines into @lines
open(ECKER,"<$ecker_file");
my @eckers=<ECKER>;
close(ECKER);

#load peak file
open(PEAKF,"<$peakfilenm");
my @peaks=<PEAKF>;
close(PEAKF);

my $totlines=$#eckers+1; #note that this is also the sequence length
my $total_overlap = 0;

# **CALCULATE VALUE FOR EACH PEAK
my $slopval=0;
my $index = 0;

PEAK: foreach my $thispeak (@peaks)
{
	chomp($thispeak);
	my @pieces = split "\t", $thispeak;
	my $pstart=$pieces[3];
	my $pend=$pieces[4];
	SEARCHLOOP: for (my $j= $index; $j <= $totlines; $j++)
	{
		chomp($eckers[$j]);
		my @pts_e = split "\t", $eckers[$j];
		my $position =$pts_e[3];
		if ($pieces[0] eq $pts_e[0])
		{	
			if( ($position >= $pstart) and ($position <= $pend) )
			{
				$nums{$pts_e[0]}->{$pstart}++;
				$index = $j;
			}
			elsif ($position > $pend )
			{
				last SEARCHLOOP;
			}
		}			
		elsif ($pieces[0] lt $pts_e[0])
		{
			last SEARCHLOOP;
		}
	}
}	

foreach my $chr (sort keys %nums)
{
	foreach my $pos (sort {$a <=> $b} keys %{$nums{$chr}})
	{
		print "$chr\t$pos\t$nums{$chr}->{$pos}\n";
	}
}

print "\a";
exit;

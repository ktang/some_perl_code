#!/usr/bin/perl -w
# this is version 5, the other 4 version before
#does't use the score in the Ecker data.
use strict;

#START HERE MAIN PROGRAM

my $usage = "$0 <Ecker_data_sorted> <peaks_gff_sorted>";
die $usage unless (@ARGV == 2);


my $ecker_file = $ARGV[0];
my $peakfilenm=$ARGV[1];
my %nums;



#load gff lines into @lines
open(ECKER,"<$ecker_file")
	or die "cannot open file $ecker_file:$!";

#load peak file
open(PEAKF,"<$peakfilenm")
	or die "cannot open file $peakfilenm:$!";

my @eckers=<ECKER>;
close(ECKER);
my @peaks=<PEAKF>;
close(PEAKF);

my $c_num = $#eckers + 1; #note that this is also the sequence length
my $p_num = $#peaks + 1;
my $total_overlap = 0;

# **CALCULATE VALUE FOR EACH PEAK
my $slopval=0;
my ($p,$c) = (0 , 0);

LOOP: while ( ($p < $p_num) and ($c < $c_num))
{
	my $thispeak = $peaks[$p];
	my $thisC = $eckers[$c];
	chomp $thispeak;
	chomp $thisC;
	my @pts_p = split "\t", $thispeak;
	my @pts_c = split "\t", $thisC;
	if( $pts_p[0] eq $pts_c[0])
	{
		if( ($pts_c[3] >= $pts_p[3] )and($pts_c[3] <= $pts_p[4] ) )
		{
			$nums{$pts_p[0]}->{$pts_p[3]} += $pts_c[6] ;
			$c++;
			$total_overlap++;
			next LOOP;		
		}
		elsif ($pts_c[3] < $pts_p[3] ) 
		{
			$c++;
			next LOOP;
		}
		elsif($pts_c[3] > $pts_p[4] )
		{
			$p++;
			next LOOP;
		}
	}

	elsif ( $pts_p[0] lt $pts_c[0] ) 
	{
		$p++;
		next LOOP;
	}
	
	elsif ( $pts_p[0] gt $pts_c[0] )
	{
		$c++;
		next LOOP;
	}

}

foreach my $chr (sort keys %nums)
{
	foreach my $pos (sort {$a <=> $b} keys %{$nums{$chr}})
	{
		my $val = $nums{$chr}->{$pos}/100;
		print "$chr\t$pos\t$val\n";
	}
}

print STDERR  "\a";
exit;

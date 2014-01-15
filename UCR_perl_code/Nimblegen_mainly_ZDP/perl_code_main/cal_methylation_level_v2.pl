#!/usr/bin/perl

# this version use the input CG_count_wsize to calculate 
# the % methylated CG.

use warnings;
use strict;

my $usage = "<$0> <Ecker_data> <CG_count_data> <output>";

die $usage unless(@ARGV == 3);

open (ECKER, "<$ARGV[0]")
	or die "cannot open Ecker_data $ARGV[0]:$!";

open (CG, "<$ARGV[1]")
	or die "cannot open input file:$!";
	
open (OUT ,">$ARGV[2]")
	or die "cannot open output file:$!";

my $line;
my @pts;
my %nums;
my @ecker;
my @CG_count;
my $index;

@ecker = <ECKER>;
close (ECKER);

@CG_count = <CG>;
close (CG);

my $ecker_num = $#ecker + 1;
my $CG_count_num = $#CG_count + 1;

my ($e, $cg) = (0, 0);

LOOP : while ( ($e < $ecker_num) and ($cg < $CG_count_num))
{
	my $thisecker = $ecker[$e];
	my $thisCG = $CG_count[$cg];
	chomp $thisecker;
	my @pts_e = split "\t",$thisecker;
	chomp $thisCG;
	my @pts_cg = split "\t", $thisCG;
	my $pos = $pts_e[3];

	if ($pts_e[0] eq $pts_cg[0])
	{
		if ( ($pos >= $pts_cg[3]) and ($pos <= $pts_cg[4]))
	       	{
			if(($pts_e[5] eq "CG") and ($pts_e[6] > 0))
			{
				$nums{$pts_cg[0]}-> {$pts_cg[3]}++; 
			}
		
			$e++;
			next LOOP;
		}
		
		elsif ($pos < $pts_cg[3] )
		{
			$e++;
			next LOOP;
		}
		elsif ( $pos > $pts_cg[4])
		{
			$cg++;
			next LOOP;
		}
	}

	elsif( $pts_cg[0] lt $pts_e[0])
	{
		$cg++;
		next LOOP;			
	}

	elsif ( $pts_cg[0] gt $pts_e[0])
	{
		$e++;
		next LOOP;
	}
}

=pod
foreach my $chr (sort keys %nums)
{
	foreach my $pos (sort {$a <=> $b} keys %{$nums{$chr}})
	{
		my $start = $pos * $wsize;
		my $end = ($pos + 1) * $wsize - 1;
		print OUT "$chr\t$start\t$end\t$nums{$chr}->{$pos}\n";
	}
}
=cut

foreach my $myline (@CG_count)
{
	my $per;
	chomp $myline;
	@pts = split "\t",$myline;
	if ( (exists $nums{$pts[0]}) and (exists $nums{$pts[0]}->{$pts[3]}))
	{$pts[6] = $nums{$pts[0]}->{$pts[3]} }
	else {$pts[6] = 0 };
	if( $pts[5] > 0)
	{
		 $per = $pts[6]/$pts[5];
	}
	elsif ($pts[6] > 0)
	{
		print "Error:$myline,$pts[6]\n";
	}
	else 
	{
		$per = 0;
	}
	$pts[7] = $per;
	my $out_line = join("\t",@pts);
	print OUT "$out_line\n";
}

close(OUT);
print STDERR"\a";
exit;

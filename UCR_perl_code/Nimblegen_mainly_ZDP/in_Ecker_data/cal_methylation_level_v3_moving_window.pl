#!/usr/bin/perl

# this version use the input CG_count_wsize to calculate 
# the % methylated CG.

use warnings;
use strict;

my $usage = "<$0> <Ecker_data> <CG_count_data> <output> <move_size>";

die $usage unless(@ARGV == 4);

open (ECKER, "<$ARGV[0]")
	or die "cannot open Ecker_data $ARGV[0]:$!";

open (CG, "<$ARGV[1]")
	or die "cannot open input file:$!";
	
open (OUT ,">$ARGV[2]")
	or die "cannot open output file:$!";

my $move_size = $ARGV[3];
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
my $thisEcker;
my $ecker_num = $#ecker + 1;
my $CG_count_num = $#CG_count + 1;

my ($e, $cg) = (0, 0);
my $former_chr = "chr0";
my ($CG_pre, $former_CG_pre,$CG_post, $former_CG_total) = (0,0,0,0) ;
my ($total,$former_total) =(0,0);
my $CG_mid = 0;
LOOP1 : while ($cg < $CG_count_num)
{
	my $thisCG = $CG_count[$cg];
	chomp $thisCG;
	my @pts_cg = split "\t", $thisCG;
		
	my $win_start = $pts_cg[3];
	my $win_end = $pts_cg[4];
	my $mid_1 = $win_start + $move_size;
	my $mid_2 = $win_end - $move_size;

# the first time of a new chromosome
   if ( $pts_cg[0] ne $former_chr )
   {
		$former_chr = $pts_cg[0];
		($CG_pre, $former_CG_pre,$CG_post, $former_CG_total) = (0,0,0,0) ;
		($total,$former_total) =(0,0);
		 $CG_mid = 0;
	LOOP1 : while ($cg < $CG_count_num)
	{
		LOOP2:	while (1)
		{
			$thisEcker = $ecker[$e];
			chomp $thisEcker;
			my @pts_e = split "\t", $thisEcker;
			my $position = $pts_e[3];			
			if ( $position < $mid_1  )
			{
				if(($pts_e[5] eq "CG") and ($pts_e[6] > 0))
				{
					$CG_pre++;
				}
				$e++;
				next LOOP2;
			}
			elsif ( $position < $mid_2 )
			{
				if(($pts_e[5] eq "CG") and ($pts_e[6] > 0))
				{
					$CG_mid++;
				}
				$e++;
				next LOOP2;
			}
			elsif ( $position <= $win_end)
			{
				if(($pts_e[5] eq "CG") and ($pts_e[6] > 0))
				{
					$CG_post++;
				}
				$e++;
				next LOOP2;
			}
			else 
			{
				$total = $CG_pre + $CG_mid + $CG_post;
				$nums{$former_chr}->{$win_start} = $total;
				last LOOP2;
			}
		}
	}
  }
	else 
	{
		$thisEcker = $ecker[$e];
		chomp $thisEcker;
		my @pts_e = split "\t", $thisEcker;
		my $position = $pts_e[3];			
		LOOP3: while(1)
		{
			if ( ($position >= $mid_2) and  ($position <= $win_end) )
			{
				if(($pts_e[5] eq "CG") and ($pts_e[6] > 0))
				{
					$CG_post++;
				}
				$e++;
				next LOOP3;
			}		
			else 
			{
				$total = $former_total - $former_CG_pre + $CG_post;
				$nums{$former_chr}->{$win_start} = $total;
				last LOOP3;
			}		
		}
	}

	$former_CG_pre = $CG_pre;
	$former_total = $total;
	$CG_pre = 0;
	$CG_mid = 0;
	$CG_post = 0;
	$cg++;
	next LOOP1;
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

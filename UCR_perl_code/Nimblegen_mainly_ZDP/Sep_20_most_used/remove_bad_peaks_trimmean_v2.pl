#! /usr/bin/perl -w

use strict;

=head
this script is used to remove the
peaks which the trim mean of orignal intensity is
negative.
=cut


sub trimmean
{
	my $arr_ref = $_[0];
	my $length = scalar (@$arr_ref);
	
	if ($length < 2)
	{
		print STDERR join "\t",@{$arr_ref};
		die "@$arr_ref";	
	}
	
	my $sum = 0;
	for(my $i = 1; $i < $length-1 ; $i++)
	{
		$sum += $arr_ref->[$i];	
	}	
	my $trimmean = $sum/($length - 2);
	return $trimmean;
}

##note : all the gff should be sgr
my $usage = "$0 <bed_file> <mutant_norm_sgr1> <mutant_norm_sgr2> <cutoff>";

die $usage unless (@ARGV == 4);

my $bed_file = $ARGV[0];

my $gff1_file = $ARGV[1];

my $gff2_file = $ARGV[2];

my $cutoff = $ARGV[3];

open (BED, "$bed_file")
	or die "cannot open $bed_file:$!";
	
open (G1, "$gff1_file")
	or die "cannot open $gff1_file:$!";
	
open (G2, "$gff2_file")
	or die "cannot open $gff2_file:$!";
	
my @beds = <BED>;
close (BED);


my @gff1 = <G1>;
close (G1);

my @gff2 = <G2>;
close (G2);

#make a flag for every peak, initiliazed as 1
my %flag;
my %scores1;
my %scores2;

print STDERR "debug";

for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my ($chr,$start, $end ) = split "\t",$this;
	$flag{$chr}->{$start} = 1;
	$scores1{$chr}->{$start} = [];
	$scores2{$chr}->{$start} = [];
}


my $b = 0; #index to go through @beds
my $g = 0; #index to go though @gff1/2

LOOP: while (($b <= $#beds) and ($g <= $#gff1))
{
	my $this_gff1 = $gff1[$g];
	my $this_gff2 = $gff2[$g];
	my $this_bed = $beds[$b];
	chomp $this_gff1;
	chomp $this_gff2;
	chomp $this_bed;
	my @pts_g1 = split "\t",$this_gff1;
	my @pts_g2 = split "\t",$this_gff2;
	my @pts_bed = split "\t", $this_bed;
	my $pos = $pts_g1[1] + 0;
	my $peak_start = $pts_bed[1] + 0;
	my $peak_end = $pts_bed[2] + 0;
	my $val1 = $pts_g1[2];
	my $val2 = $pts_g2[2];
	
	my $sgr1_pos = $pts_g1[1] + 0;
	my $sgr2_pos = $pts_g2[1] + 0;

	if (( $pts_g1[0] ne $pts_g2[0]) or ($sgr1_pos != $sgr2_pos))
	{
		print STDERR "$this_gff1\n$this_gff2\n";
		die "not the same line in sgrs";
	}
	
	if ($pts_g1[0] eq $pts_bed[0])
	{
		if( $pos < $peak_start )
		{
			$g++;
			next LOOP;	
		}
		
		elsif ( ($pos >= $peak_start) and ($pos <= $peak_end ))
		{

############################
=head1
			if( ($val1 < 0) or ($val2 < 0))
			{
				if ( $flag{$pts_bed[0]}->{$pts_bed[1]} == 1)
				{$flag{$pts_bed[0]}->{$pts_bed[1]} = -1;}
				elsif ( $flag{$pts_bed[0]}->{$pts_bed[1]} == -1)
				{}
				else {die "hash wrong";}
				$b++;
				next LOOP;
			}
=cut
############################			
			push @{$scores1{$pts_bed[0]}->{$pts_bed[1]}},$val1;
			push @{$scores2{$pts_bed[0]}->{$pts_bed[1]}},$val2;

			$g++;
			next LOOP;
		}
		
		elsif ($pos > $peak_end )
		{
			$b++;
			next LOOP;	
		}
		else 
		{
			print STDERR "$this_bed,$this_gff1,$this_gff2";
			die "same chr but position is impossible";
		}
	}
	elsif ( $pts_g1[0] gt $pts_bed[0])
	{
		$b++;
		next LOOP;	
	}
	elsif ($pts_g1[0] lt $pts_bed[0] )
	{
		$g++;
		next LOOP;	
	}
	else
	{
		print STDERR "$this_bed,$this_gff1,$this_gff2";
		die "impossible chr";
	}
}

my $num_less_four;

for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my ($chr,$start, $end ) = split "\t",$this;
	
	my @array1 = sort @{$scores1{$chr}->{$start}};
	my @array2 = sort @{$scores2{$chr}->{$start}};
	
	if( ($#array1 < 4)or ($#array2 < 4))
	{
			$num_less_four++;
	}
	else
	{
		my $mean_1 = trimmean(\@array1);
		my $mean_2 = trimmean(\@array2);
	
		if ( ($mean_1 >= $cutoff) and ($mean_2 >= $cutoff) )
		{
			print  "$this\n";	
		}
	}
}

print STDERR "num_less_four = $num_less_four\n";
exit;
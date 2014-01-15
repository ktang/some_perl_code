#! /usr/bin/perl -w

use strict;

=head
this script is used to remove the
peaks which the orignal intensity is
negative.

=cut


=head
sub isoverlappingslop
{
	my $extrabit=$_[0]; #allows for sloppiness - if this is set to zero, then there is no sloppiness
	my $acor=$_[1]-$extrabit;
	my $bcor=$_[2]+$extrabit;
	my $ccor=$_[3]-$extrabit;
	my $dcor=$_[4]+$extrabit;
	
	if ($acor>=$ccor && $acor<=$dcor)
	{return 1;}
	if ($bcor>=$ccor && $bcor<=$dcor)
	{return 1;}
	if ($ccor>=$acor && $ccor<=$bcor)
	{return 1;}
	if ($dcor>=$acor && $dcor<=$bcor)
	{return 1;}
	return 0;
}
=cut
##note : all the gff should be sgr
my $usage = "$0 <bed_file> <mutant_norm_sgr1> <mutant_norm_sgr2>";

die $usage unless (@ARGV == 3);

my $bed_file = $ARGV[0];

my $gff1_file = $ARGV[1];

my $gff2_file = $ARGV[2];

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

for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my ($chr,$start, $end ) = split "\t",$this;
	$flag{$chr}->{$start} = 1;
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

for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my ($chr,$start, $end ) = split "\t",$this;
	if ($flag{$chr}->{$start} == -1)
	{
		print STDERR "$this\n";	
	}
	elsif ($flag{$chr}->{$start} == 1)
	{
		print "$this\n";
	}
	else 
	{
		print STDERR "wrong flag";
		die;	
	}
}


#print STDERR "\a";
exit;
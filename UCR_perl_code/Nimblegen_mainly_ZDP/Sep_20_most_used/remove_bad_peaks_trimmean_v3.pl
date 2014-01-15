#! /usr/bin/perl -w

use strict;

=head
this script is used to remove the
peaks which the trim mean of orignal intensity is
negative.and the output is annoted.
=cut


sub trimmean
{
	my $arr_ref = $_[0];
	my $length = scalar (@$arr_ref);
	
	if ($length < 2)
	{
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
my $usage = "$0 <bed_file> <mutant_norm_sgr1> <mutant_norm_sgr2> <cutoff> <label_or_hypothesis>";

die $usage unless (@ARGV == 5);

my $sybol = $ARGV[4];
if ($sybol ne 'label' and $sybol ne 'hypothesis')
{
	die $usage;	
}

my $bed_file = $ARGV[0];

my $gff1_file = $ARGV[1];

my $gff2_file = $ARGV[2];

my $cutoff = $ARGV[3];

my $annoted_file1 = "/Users/kaitang/Desktop/Nimblegen/two_list_for_ZDP_wrong_analysis/send/labeled_ZDP_peaks_with_anno.txt";

my $annoted_file2 = "/Users/kaitang/Desktop/Nimblegen/two_list_for_ZDP_wrong_analysis/send/hypothesis_ZDP_peaks_with_anno.txt";

open (LA, "$annoted_file1")
	or die "cannot open $annoted_file1";
	
open (HY, "$annoted_file2")
	or die "cannot open $annoted_file2";
	
open (BED, "$bed_file")
	or die "cannot open $bed_file:$!";
	
open (G1, "$gff1_file")
	or die "cannot open $gff1_file:$!";
	
open (G2, "$gff2_file")
	or die "cannot open $gff2_file:$!";
	
my @las = <LA>;
close(LA);
my @hys = <HY>;
close (HY);

my @results;

if ($sybol eq "label")
{
	@results = @las;
}

elsif ($sybol eq "hypothesis")
{
	@results = @hys;
}
	
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

my %outs;

for (my $i = 0; $i <= $#beds; $i++)
{
	my $this = $beds[$i];
	chomp $this;
	my ($chr,$start, $end ) = split "\t",$this;
	
	my @array1 = sort @{$scores1{$chr}->{$start}};
	my @array2 = sort @{$scores2{$chr}->{$start}};
	
	my $mean_1 = trimmean(\@array1);
	my $mean_2 = trimmean(\@array2);
	
	if ( ($mean_1 >= $cutoff) and ($mean_2 >= $cutoff) )
	{
		#print  "$this\n";
		my $mean = 0.5*($mean_1+$mean_2);
		$outs{$chr}->{$start} = $mean;
			
	}
	else
	{
		#print STDERR "$this\n";
	}
}

for (my $i = 0; $i <= $#results; $i++)
{
	my $this = $results[$i];
	chomp $this;
	my @pts = split "\t", $this;
	if (exists $outs{$pts[0]}->{$pts[1]})
	{
		$pts[3] = $outs{$pts[0]}->{$pts[1]};
		my $out_line = join "\t",@pts;
		print $out_line,"\n";
	}	
}

=head2
for (my $i = 0; $i <= $#las; $i++)
{
	my $this = $las[$i];
	chomp $this;
	my @pts = split "\t", $this;
	if (exists $outs{$pts[0]}->{$pts[1]})
	{
		$pts[3] = $outs{$pts[0]}->{$pts[1]};
		my $out_line = join "\t",@pts;
		print $out_line,"\n";
	}	
}

for (my $i = 0; $i <= $#hys; $i++)
{
	my $this = $hys[$i];
	chomp $this;
	my @pts = split "\t", $this;
	if (exists $outs{$pts[0]}->{$pts[1]})
	{
		$pts[3] = $outs{$pts[0]}->{$pts[1]};
		my $out_line = join "\t",@pts;
		print $out_line,"\n";	
	}	
}
=cut

#print STDERR "\a";
exit;
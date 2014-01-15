#!/usr/bin/perl -w

=head1
use file sorted by signal as input
parameter num_of_level, num_per_level
=cut

use strict;

my $usage = "$0 <input_gff_sorted_by_signal> <num_of_level> <num_per_level>";

die $usage unless (@ARGV == 3);

my $level_num = $ARGV[1];

my $num_per_level = $ARGV[2];

open (GFF,"$ARGV[0]")
	or die "cannot opne $ARGV[0]:$!";
my @gffs = <GFF>;
close(GFF);

my $total = $#gffs + 1;

my $interval = int ($total/($level_num - 1));

for(my $i = 0; $i < $num_per_level; $i ++)
{
	my $this= $gffs[$i];
	chomp $this;
	my @pts = split"\t", $this;
	$pts[8] = $i + 1;
	my $out = join "\t", @pts;
	print "$out\n";	
}

for (my $i = 1; $i <= $level_num - 2; $i++)
{
	my $start = $i * $interval;
		for (my $j = 0; $j < 100; $j ++)
		{
			my $this= $gffs[$start + $j];
			chomp $this;
			my @pts = split"\t", $this;
			$pts[8] = $start + $j + 1;
			my $out = join "\t", @pts;
			print "$out\n";	
		}
}

for(my $i = $total - $num_per_level; $i < $total; $i ++)
{
	my $this= $gffs[$i];
	chomp $this;
	my @pts = split"\t", $this;
	$pts[8] = $i + 1;
	my $out = join "\t", @pts;
	print "$out\n";		
}

exit;
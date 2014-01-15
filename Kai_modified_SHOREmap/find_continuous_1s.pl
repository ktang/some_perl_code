#!/usr/bin/perl -w

use strict;

my $usage = "$0 <input> <cutoff> <min_continuous_number>";

die $usage unless(@ARGV == 3);

my $input = $ARGV[0];
my $cutoff = $ARGV[1];
my $minlen = $ARGV[2];



open(IN,"<$input");
my @lines=<IN>;
close(IN);

my $inpiece=0;
my $startpoint=0;
my $endpoint=0;
my $lenpeak=0;
my $piece_region_name="";

LOOPER: foreach my $thisone (@lines)
{
	chomp($thisone);
	my @pieces=split "\t", $thisone;
	if ($inpiece and $piece_region_name ne $pieces[0])
	{	#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print join("\t", ($piece_region_name, $startpoint, $endpoint)), "\n";
		}
		$lenpeak=0;
		#note it is important that there is no "next LOOPER" statement here because this point could be the start of a new
		# hit
	}
	
	if ($inpiece==0 and $pieces[10] > $cutoff)
	{
		next LOOPER;
	}
	if ($inpiece==0 and $pieces[10] <= $cutoff)
	{
			
		$startpoint=$pieces[1];
		$endpoint=$pieces[1]; #normally, this would not be endpoint of run but in weird case where we want single point runs, this allows it
		$inpiece=1;
		$lenpeak=1;
		$piece_region_name=$pieces[0];
		next LOOPER;
	}
	
	if ($inpiece and $pieces[10] <= $cutoff)
	{
		$endpoint=$pieces[1];
		$lenpeak++;
		next LOOPER;
	}
	if ($inpiece and $pieces[10] >  $cutoff)
	{
		#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print join("\t", ($piece_region_name, $startpoint, $endpoint)) , "\n";
		}
		$lenpeak=0;
		next LOOPER;
	}
}

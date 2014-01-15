#!/usr/bin/perl -w

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

#short for bed, long for gff
my $usage = "$0 <bed_file> <gff_file> <gap>";
die $usage unless(@ARGV == 3);

my $fname1=$ARGV[0]; #input filenam
my $fname2=$ARGV[1]; #second filename
my $slopval=$ARGV[2] + 0; #extra bit to add to each piece to allow near overlaps

my @sps;
my @lps;
my $line;
my ($send,$sstart,$lend,$lstart);
my ($long_num,$short_num,$overlap) =(0,0,0);

print "$long_num,$short_num,$overlap\n";

open (SHORT,"<$fname1")
	or die "cannot open file $fname1:$!";
open (LONG,"<$fname2")
	or die "cannot opne file $fname2:$!";

my @short=<SHORT>;
close(SHORT);

$short_num = $#short +1;

while ($line = <LONG>)
{		
		$long_num++;
		chomp($line);
		@lps = split "\t", $line;
	LOOKLOOP: for ($j=0; $j <= $#short; $j++)
	{
		@sps = split "\t", $short[$j];
		$sstart=$sps[1];
		$send=$sps[2];
		$lstart = $lps[3];
		$lend = $lps[4];
		if (($lps[0] eq $sps[0]) and isoverlappingslop($slopval, $lstart,$lend,$sstart,$send))
		{
			$overlap++;
			last LOOKLOOP;
		}
	}
}
close(LONG);

my $gaps = 2* $slopval;

print "long file contains $long_num lines \n short files contains $short_num lines\n overlap within $gaps bp is $overlap\n\a\n";
exit;

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
my $usage = "$0 <anno_file> <compared_bed_file> <gap>";
die $usage unless(@ARGV == 3);

my $annotated_file = $ARGV[0]; #input annotated file
my $col_file=$ARGV[1]; #other compare file, add as a col if has overlap with $annotated_file
my $slopval=$ARGV[2] + 0; #extra bit to add to each piece to allow near overlaps

my @sps;
my @lps;
my $line;
my ($send,$sstart,$lend,$lstart);
my ($long_num,$short_num,$overlap) =(0,0,0);

print STDERR "$long_num,$short_num,$overlap\n";

open (ANNO,$annotated_file)
	or die "cannot open file $annotated_file:$!";
open (COL, $col_file)
	or die "cannot opne file $col_file:$!";

my @annos=<ANNO>;
close(ANNO);

$short_num = $#annos +1;

#long is col_file
#short is annotated_file

my %hashs;

while ($line = <COL>)
{		
		$long_num++;
		chomp($line);
		@lps = split "\t", $line;
	LOOKLOOP: for ($j=0; $j <= $#annos; $j++)
	{
		@sps = split "\t", $annos[$j];
		$sstart=$sps[1];
		$send=$sps[2];
		$lstart = $lps[1];
		$lend = $lps[2];
		if (($lps[0] eq $sps[0]) and isoverlappingslop($slopval, $lstart,$lend,$sstart,$send))
		{
			$overlap++;
			$hashs{$sps[0]}->{$sps[1]} = $lps[0].":".$lps[1]."-".$lps[2];
			last LOOKLOOP;
		}
	}
}
close (COL);

for (my $i = 0; $i <= $#annos; $i++){
	my $this = $annos[$i];
	chomp $this;
	my @pts = split "\t",$this;
	my $col = "NONE";
	if (defined $hashs{$pts[0]}->{$pts[1]}){
		$col = $hashs{$pts[0]}->{$pts[1]};
	}
	print join("\t",(@pts,$col)),"\n";
}

my $gaps = 2* $slopval;

print STDERR "long file contains $long_num lines \n short files contains $short_num lines\n overlap within $gaps bp is $overlap\n\a\n";
exit;

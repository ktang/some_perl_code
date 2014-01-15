#!/usr/bin/perl -w
# user input an int, count the original pileup file
# to find how many positions depth is >= x

my $usage ="$0 <input_pileup> <int_cutoff>";

die $usage unless (@ARGV == 2);

open (IN,"<$ARGV[0]")
	or die "cannot open input file $ARGV[0]:$!";
	
my $line;
my $count = 0;
my @parts;

while ($line =  <IN>)
{
	chomp $line;
	@parts = split /\t/,$line;
	if ($parts[3] >= $ARGV[1])
	{$count++}
}


print $count,"\n";

  print "\a";
  exit;

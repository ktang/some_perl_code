#!/usr/bin/perl -w

# remove lines which all the depth are 0, and also
# add a column total

use strict;

my $usage = "$0 <cut_table> <output>";

die $usage unless(@ARGV == 2);

open (IN,"<$ARGV[0]") or die "Cannot open the input file $ARGV[0]:$!";

open (OUT,">$ARGV[1]") or die "Cannot open the output file $ARGV[1]:$!";

my $line;
my @parts;

$line = <IN>;
chomp $line;

print OUT "$line\ttotal\n";

while ($line = <IN>)
{
	chomp $line;
	@parts = split /\t/, $line;
	my $total = $parts[4]+$parts[5]+$parts[6]+$parts[7]+$parts[8]+$parts[9]+$parts[10]
	+$parts[11]+$parts[12]+$parts[13]+$parts[14]+$parts[15]+
	$parts[16]+$parts[17]+$parts[18]+$parts[19];
	
	if ( $total == 0 )
	{;}
	else 
	{
		print OUT $line,"\t$total\n";
	}
}


close IN;
close OUT;
exit;
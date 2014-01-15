#!/usr/bin/perl -w

use strict;

=head
This program calculate the difference between the 
two files, the format is chr position number_of_mC 
=cut

##begin of main

my $usage = "$0 <input_WT> <input_mutant>";

die $usage unless (@ARGV == 2);

open (WT, "<$ARGV[0]")
	or die "cannot open the input file $ARGV[0]:$!";
	
open (MU, "<$ARGV[1]")
	or die "cannot open the file $ARGV[1]:$!";

	
my %wts;
my $line;
my @pts;
my $diff;

while ($line = <WT>)
{
	chomp $line;
	@pts = split "\t" , $line;
	$wts{$pts[0]}->{$pts[1]} = $pts[2];	
}
close(WT);

while ($line = <MU>)
{
	chomp $line;
	@pts = split "\t", $line;
	if (exists $wts{$pts[0]} and (exists $wts{$pts[0]}->{$pts[1]} ))
	{
		
		$diff = $pts[2] - $wts{$pts[0]}->{$pts[1]} ;
		if ($diff > 0)
		{
			$pts[2] = $diff;
			$line = join ("\t", @pts);
			print   "$line\n";	
		}
			
	}	
	else 
	{
		print  "$line\n";	
	}
}

close (MU);
print STDERR "\a";
exit;

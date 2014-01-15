#!/usr/bin/perl -w

use strict;




=head
This program calculate the difference between the 
rdd methylation file and the Col-0.
=cut

##begin of main

#input

my $usage = "$0 <input> <output>";

die $usage unless (@ARGV == 2);

open (IN, "<$ARGV[0]")
	or die "cannot open the input file $ARGV[0]:$!";
	
open (WT, "<./tair9_mc_col0_wsize_500")
	or die "cannot open the file tair9_mc_col0_wsize_500:$!";

open (OUT,">$ARGV[1]")
	or die "$!";
	
my %wts;
my $line;
my @pts;
my $diff;

while ($line = <WT>)
{
	chomp $line;
	@pts = split "\t" , $line;
	$wts{$pts[0]}->{$pts[1]} = $pts[3];	
}
close(WT);

#my $i =0;

while ($line = <IN>)
{
	
#	$i ++;
	chomp $line;
	@pts = split "\t", $line;
	if (exists $wts{$pts[0]} and (exists $wts{$pts[0]}->{$pts[1]} ))
	{
		
		
						
			$diff = $pts[3] - $wts{$pts[0]}->{$pts[1]} ;
			
#if($i <30)
#			{
#				print "$pts[3] \t $wts{$pts[0]}->{$pts[1]}\n"
#			}
			
			if ($diff > 0)
			{
				$pts[3] = $diff;
				$line = join ("\t", @pts);
				print  OUT "$line\n";	
			}
			
	}	
	else 
	{
		print  OUT "$line\n";	
	}
}

close (IN);
close (OUT);
print "\n\a\a\n";
exit;

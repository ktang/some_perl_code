#!/usr/bin/perl -w

# extract pos_C24_freq.txt with at least a depth >= 5

use strict;


my $usage = "$0 <input> <output>";

die $usage unless(@ARGV == 2);

my $line;
my @parts;

open (IN, "<$ARGV[0]")
	or die "Cannot open $ARGV[0]:$!";
	
open (OUT, ">$ARGV[1]")
  or die "Cannot open $ARGV[2]: $!";
  
 print OUT "gene\tposition\tCol0_base\tC24_base\tdepth_Col0_flowcell54_lane1\tdepth_C24_flowcell54_lane1\tfreq_C24_flowcell54_lane1\tdepth_Col0_flowcell54_lane2\tdepth_C24_flowcell54_lane2\tfreq_C24_flowcell54_lane2\tdepth_Col0_flowcell54_lane3\tdepth_C24_flowcell54_lane3\tfreq_C24_flowcell54_lane3\tdepth_Col0_flowcell54_lane5\tdepth_C24_flowcell54_lane5\tfreq_C24_flowcell54_lane5\tdepth_Col0_flowcell54_lane6\tdepth_C24_flowcell54_lane6\tfreq_C24_flowcell54_lane6\tdepth_Col0_flowcell54_lane7\tdepth_C24_flowcell54_lane7\tfreq_C24_flowcell54_lane7\tdepth_Col0_flowcell54_lane8\tdepth_C24_flowcell54_lane8\tfreq_C24_flowcell54_lane8\tdepth_Col0_flowcell55_lane1\tdepth_C24_flowcell55_lane1\tfreq_C24_flowcell55_lane1\n";

while ($line = <IN>)
{
	chomp($line);
	@parts = split /\t/,$line;
	if( $parts[4] >= 5 ||$parts[5] >= 5 ||$parts[7] >= 5 ||$parts[8] >= 5 ||$parts[10] >= 5 ||$parts[11] >= 5 ||$parts[13] >= 5 ||$parts[14] >= 5 ||$parts[16] >= 5 ||$parts[17] >= 5 ||$parts[19] >= 5 ||
	$parts[20] >= 5 ||$parts[22] >= 5 ||$parts[23] >= 5 ||$parts[25] >= 5 ||$parts[26] >= 5 )
	{print OUT $line,"\n";}
	
}


close (IN);
close(OUT);
 
 print "\a";
 exit;
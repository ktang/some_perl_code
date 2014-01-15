#!/usr/bin/perl


=head2

PROGRAM: threshpeakv03.pl
GOAL: this is designed for analysis of tiling array data
this program takes a file, a threshold for hits and a minimum peak width and finds all the peaks in the data

INPUT: from command line, three pieces of data
threshold min_len filename
threshold is threshold for value in data
min_len is minimum number of consecutive points that must be above threshold to count as a peak
filename is the name of the file: MUST BE GFF FORMAT

OUTPUT: to STDOUT
gff format file

NOTE FOR  USE WITH ENCODE PIECES
(1) to maximize robustness and avoid accidentally peaks crossing boundaries of ENCODE regions,
"encode region" coding should be used for lines - like ENr323 instead of chr17 (chr17 may be wrong)

TESTING STATUS
not tested yet at all...

OLDER PROGRAM NOTES
#command line is <threshold> <min len> <filename>

#this one is designed so that if the region name changes, then there is a break in the peak-finding
#this is critical because imagine that there is a run of above threshold points at the end of one region and the start of another - 
#this program would put those together incorrectly...

VERSION NOTES
threshpeakv03.pl
(1) eliminates useless stuff with firstpiece variable
(2) makes slightly more robust - will handle single point runs if necessary.

=cut


$thresh=$ARGV[0];
$minlen=$ARGV[1];
$filenm=$ARGV[2];

open(GFF,"<$filenm");
@lines=<GFF>;
close(GFF);

$inpiece=0;
$startpoint=0;
$endpoint=0;
$lenpeak=0;
$piece_region_name="";

LOOPER: foreach $thisone (@lines)
{
	chomp($thisone);
	@pieces=split "\t", $thisone;
	if ($inpiece and $piece_region_name ne $pieces[0])
	{	#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print "$piece_region_name\tthreshpeakv02.pl\t$pieces[2]\t$startpoint\t$endpoint\t$lenpeak\t+\t.\t$piece_last_field\n";
		}
		$lenpeak=0;
		#note it is important that there is no "next LOOPER" statement here because this point could be the start of a new
		# hit
	}
	
	if ($inpiece==0 and $pieces[5]<$thresh)
	{
		next LOOPER;
	}
	if ($inpiece==0 and $pieces[5]>=$thresh)
	{
			
		$startpoint=$pieces[3];
		$endpoint=$pieces[4]; #normally, this would not be endpoint of run but in weird case where we want single point runs, this allows it
		$inpiece=1;
		$lenpeak=1;
		$piece_region_name=$pieces[0];
		$piece_last_field=$pieces[-1] . ";thresh=$thresh;minwidth=$minlen;filename=$filenm";
		next LOOPER;
	}
	
	if ($inpiece and $pieces[5]>=$thresh)
	{
		$endpoint=$pieces[4];
		$lenpeak++;
		next LOOPER;
	}
	if ($inpiece and $pieces[5]<$thresh)
	{
		#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print "$piece_region_name\tthreshpeakv03.pl\t$pieces[2]\t$startpoint\t$endpoint\t$lenpeak\t+\t.\t$piece_last_field\n";
		}
		$lenpeak=0;
		next LOOPER;
	}
}

			
	
#!/usr/bin/perl


=head2

PROGRAM: threshpeakLIMITv01.pl
GOAL: this is designed for analysis of tiling array data; this one is appropriate for CHR-based data.
this program takes a file, a threshold for hits and a minimum peak width and finds all the peaks in the data.
Note that this version has important changes from the previous version (threshpeakv03.pl).
In particular, there is max limit for the size of a gap in the data that is allowed for a peak to cross.

INPUT: from command line, three pieces of data
threshold min_len filename
threshold is threshold for value in data
min_len is minimum number of consecutive points that must be above threshold to count as a peak
filename is the name of the file: MUST BE GFF FORMAT

OUTPUT: to STDOUT
gff format file

NOTE FOR  USE WITH CHR pieces
(1) this one is specially designed to allow the use of chr coding in which there could be two widely separated
regions on a single chromosome. This one will not combine pieces across a large gap (large gap is defined by the
maxgaplength constant in the program code; currently it is 50,000 nts or 50 kb).

TESTING STATUS
previous version threshpeakv03.pl was heavily tested
however! this one REALLY NEEDS TESTING!!

OLDER PROGRAM NOTES
#command line is <threshold> <min len> <filename>

#this one is designed so that if the region name changes, then there is a break in the peak-finding
#this is critical because imagine that there is a run of above threshold points at the end of one region and the start of another - 
#this program would put those together incorrectly...

VERSION NOTES
threshpeakLIMITv01.pl
(1) uses a maximum gap length; this is good for files in which there can be multiple regions on a chromosome
threshpeakv03.pl
(1) eliminates useless stuff with firstpiece variable
(2) makes slightly more robust - will handle single point runs if necessary.

=cut

#CONSTANTS
$maxgaplength=50000; #this is the maximum length (in nucleotides) of a masked element of in the array that a peak can cross
#END OF CONSTANTS

########## add here
	
$usage ="$0 <cutoff> <number_of_continuous> <file>"

die $usage unless (@ARGV == 3)

###  end

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
	
	if ($inpiece and ($piece_region_name ne $pieces[0]))
	{	#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print "$piece_region_name\t.\t$pieces[2]\t$startpoint\t$endpoint\t$lenpeak\t+\t.\t$piece_last_field\n";
		}
		$lenpeak=0;
		#note it is important that there is no "next LOOPER" statement here because this point could be the start of a new
		# hit
	}
	if ($inpiece and (($pieces[5]<$thresh) or (($pieces[3]-$lastendpoint)>$maxgaplength))  )
	{
		#end of piece
		$inpiece=0;
		if ($lenpeak>=$minlen)
		{
			print "$piece_region_name\t.\t$pieces[2]\t$startpoint\t$endpoint\t$lenpeak\t+\t.\t$piece_last_field\n";
		}
		$lenpeak=0;
		next LOOPER;
	}
	if (($inpiece==0) and ($pieces[5]<$thresh))
	{
		next LOOPER;
	}
	if (($inpiece==0) and ($pieces[5]>=$thresh))
	{
			
		$startpoint=$pieces[3];
		$endpoint=$pieces[4]; #normally, this would not be endpoint of run but in weird case where we want single point runs, this allows it
		$lastendpoint=$endpoint;
		$inpiece=1;
		$lenpeak=1;
		$piece_region_name=$pieces[0];
		$piece_last_field=$pieces[-1] . ";thresh=$thresh;minwidth=$minlen;filename=$filenm";
		next LOOPER;
	}
	
	if ($inpiece and ($pieces[5]>=$thresh))
	{
		$endpoint=$pieces[4];
		$lastendpoint=$endpoint;
		$lenpeak++;
		next LOOPER;
	}
	
}

print STDERR "\a";			
exit;
#!/usr/bin/perl

=head2

PROGRAM: ENCODErankmax4v02.pl
GOAL: this program produces ranks for each ENCODE peak. 
IMPORTANT: for this one, the ranks are based on the "maxfour" value for each peak
(The alternate version is ENCODErankv01.pl, which calculates a maximum z-score.

INPUT: 2 files on command line
(1) file of peaks - WITH SAME ENCODING AS gff file
(2) SORTED gff file. File should be sorted by sortgfffilev02.sh BEFORE using

BIG CHANGE for v02 version:
(1) this one outputs the coordinates of the maxfour location instead of the coords of the original peak
this one writes the coordinates of the original peak into the last field of the gff
(2) the old version (v01) just used the original coordinates...

MINOR v02 CHANGE
(1) I killed a bunch of the old code that was commented out in this one
(2) this old code is still available in v01 version etc.

OUTPUT: to STDOUT

NOTES:
(1) this is sorta slow; could be faster...
(2) I made this one output to three decimal places
(3) I also fixed the maxfour stuff in this one... it was not doing it quite right. 
(4) this one also outputs stuff if commands are not listed...


SOME NOTES ON WHY I DID THIS
So I think that this one will be more robust for estimating fold enrichment... but maybe not...

NOTE ON TESTING:
essentially untested at this point; does run.

EXAMPLE
ENCODErankmax4v02.pl wga_spikeins_L2.gff ENCODE_spikeins_qmed.gff > ENCODE_spikeins_qmed_RANK.gff

=cut



#sub getratio takes a single string that is a gff line and extracts the ratio
sub getratio
{
	my $daline=$_[0];
	my @dpieces=();
	
	chomp($daline);
        	@dpieces=split "\t", $daline;
        	return ($dpieces[5]);
}




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
	
	 
	



	
#START HERE MAIN PROGRAM

#debug this is the debug flag! if it is 1, then debug output
$debug=0;

$peakfilenm=$ARGV[0];
$newgff=$ARGV[1];

if ($ARGV[1] eq "")
{print "The following is an example of usage. Note that peaks file comes FIRST and that gff input file must be sorted!!\n";
print 'ENCODErankmax4v02.pl wga_spikeins_L2.gff ENCODE_spikeins_qmed.gff > ENCODE_spikeins_qmed_RANK.gff' . "\n";}



#load gff lines into @lines
open(GFF,"<$newgff");
@lines=<GFF>;
close(GFF);

#load peak file
open(PEAKF,"<$peakfilenm");
@peaks=<PEAKF>;
close(PEAKF);



$totlines=$#lines+1; #note that this is also the sequence length



# **CALCULATE VALUE FOR EACH PEAK
$slopval=0;
foreach $thispeak (@peaks)
{


#SECTION peak-1: grab all lines from gff that are in the peak
#note - there is a chance that this grabs an extra oligo on each end...
	@peakpieces=(); #reset @peakpieces
	$foundflag=0; #this flag indicates whether it was found or not
	chomp($thispeak);
	@pieces = split "\t", $thispeak;
	$fstart=$pieces[3];
	$fend=$pieces[4];
	SEARCHLOOP: for ($j=0; $j < $totlines; $j++)
	{
		#chomp($lines[$j]);
		@spieces = split "\t", $lines[$j];
		$sstart=$spieces[3];
		$send=$spieces[4];
		#print "fstart is $fstart fend is $fend sstart is $sstart send is $send\n";
		
		
		if (($pieces[0] eq $spieces[0]) and isoverlappingslop($slopval, $fstart,$fend,$sstart,$send))
		{
			push @peakpieces, $lines[$j] ;
			$foundflag=1; #yes, I have found some!

		}
		#now, because the file is sorted, we want to escape as soon
		# as we exit 
		#this will save a huge number of searches... but maybe not so important
		#I could kill this next section if I get worried
		if ((!isoverlappingslop($slopval, $fstart,$fend,$sstart,$send)) and $foundflag)
		{last SEARCHLOOP;}
		
		
	}

#debug
if ($debug) {print "peakpieces:\n @peakpieces\n";};






#SECTION peak-2: get the value for this peak
#this now uses the maxfour value for things
#this is the section where we change so that we do the maxn stuff
#assume that there is a constant $RUNLEN that determines the runlength to test
#STOPPOINTHERE

#NOT FINISHED NOT FINISHED NOT FINISHED!!!!


	#$RUNLEN is the constant for run length
	$totpieces=$#peakpieces+1;
	$thismaxpval=-100;
	$maxav=-100;
	@ts_pieces = split "\t", $peakpieces[0];
	$maxav_start=$ts_pieces[3];
	@te_pieces = split "\t", $peakpieces[3];
	$maxav_end=$te_pieces[4];

	for ($i=0;$i<=($#peakpieces-3);$i++)
	{
		$thisav=(getratio($peakpieces[$i])+getratio($peakpieces[$i+1])+getratio($peakpieces[$i+2])+getratio($peakpieces[$i+3]))/4;
		#calculate start and end of this piece
		@ts_pieces = split "\t", $peakpieces[$i];
		$m4start=$ts_pieces[3];
		@te_pieces = split "\t", $peakpieces[$i+3];
		$m4end=$te_pieces[4];
		if (abs($thisav)<0.001)
		{$thisav=0;}
		#print "thisav is $thisav\n";
		#following line with ALLVAL added in v02 for output of every value
		print ALLVAL ($thisav . "\n");
		if ($thisav>$maxav)
		{ 	#now, readjust start and end position and actual value
			$maxav=$thisav;
		  	$maxav_start=$m4start;
		  	$maxav_end=$m4end;
		}
		$thispval=$thisav; #these are bad and incorrect names, but I am keeping them
		$thismaxpval=$maxav; #these are bad and incorrect, but I am keeping them


		#debug
		if ($debug)
		{			
		print "pointlevel, wlength, percentileval, thispval: $pointlevel, $wlength, $percentileval, $thispval\n";
		}		
		
	} #end of loop here

#SECTION peak-3: do the output for this peak
	#I've modified this to make it easier to use in signalmap
	#note that for v02 version, output coords of maxfour best signal instead of just use old
	# peak coords.
	# BUT: I retain the old coords in the 8th field
	@duppieces=@pieces;
	$duppieces[1]="ENCODErankmax4v02";
	$duppieces[2]="ENCODErankmax4v02_" . $peakfilenm; #change third field
	$duppieces[8]=$pieces[8] .';' . "origcoords=$pieces[3],$pieces[4]"; #add orig coords to last field
	$duppieces[5]=(int($thismaxpval*1000))/1000; #truncate to three decimal places
	$duppieces[3]=$maxav_start;
	$duppieces[4]=$maxav_end;
	$outputstr=join("\t",@duppieces);
	print "$outputstr\n";
	#OLD way - another OPTION
	# print $thispeak . ';ENCODErankv01max4' ."=$thismaxpval\n";
}



	 
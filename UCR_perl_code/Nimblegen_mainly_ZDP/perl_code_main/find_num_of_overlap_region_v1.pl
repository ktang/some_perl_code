#!/usr/bin/perl -w

use strict;


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

my $usage = "$0 <Ecker_data> <peaks_gff>";
die $usage unless (@ARGV == 2);

#$debug=0;

my $ecker_file = $ARGV[0];
my $peakfilenm=$ARGV[1];


#load gff lines into @lines
open(ECKER,"<$ecker_file");
my @eckers=<ECKER>;
close(ECKER);

#load peak file
open(PEAKF,"<$peakfilenm");
my @peaks=<PEAKF>;
close(PEAKF);


my $totlines=$#eckers+1; #note that this is also the sequence length
my $total_overlap = 0;

# **CALCULATE VALUE FOR EACH PEAK
my $slopval=0;


foreach my $thispeak (@peaks)
{
#SECTION peak-1: grab all lines from gff that are in the peak
#note - there is a chance that this grabs an extra oligo on each end...
	my @peakpieces=(); #reset @peakpieces
	my $foundflag=0; #this flag indicates whether it was found or not
	chomp($thispeak);
	my @pieces = split "\t", $thispeak;
	my $pstart=$pieces[3];
	my $pend=$pieces[4];
	SEARCHLOOP: for (my $j=0; $j < $totlines; $j++)
	{
		#chomp($lines[$j]);
		my @spieces = split "\t", $eckers[$j];
		my $estart=$spieces[1];
		my $eend=$spieces[2];
		#print "fstart is $fstart fend is $fend sstart is $sstart send is $send\n";
		
		if (($pieces[0] eq $spieces[0]) and isoverlappingslop($slopval, $pstart,$pend,$estart,$eend))
		{
                        my @peakpieces;
			push  @peakpieces, $eckers[$j] ;
			$total_overlap ++;
			$foundflag=1; #yes, I have found some!

		}
		#now, because the file is sorted, we want to escape as soon
		# as we exit 
		#this will save a huge number of searches... but maybe not so important
		#I could kill this next section if I get worried
		if ((!isoverlappingslop($slopval, $pstart,$pend,$estart,$eend)) and $foundflag)
		{last SEARCHLOOP;}
	}
}	

print "total overlap is $total_overlap\n\n\a";

print "\a";
exit;


###add a }


#debug
#if ($debug) {print "peakpieces:\n @peakpieces\n";};






#SECTION peak-2: get the value for this peak
#this now uses the maxfour value for things
#this is the section where we change so that we do the maxn stuff
#assume that there is a constant $RUNLEN that determines the runlength to test
#STOPPOINTHERE

#NOT FINISHED NOT FINISHED NOT FINISHED!!!!
=pod

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
=cut
#}



	 

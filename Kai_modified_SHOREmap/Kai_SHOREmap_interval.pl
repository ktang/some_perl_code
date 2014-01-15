#!/usr/bin/perl

# --------------------------------------------------------------------
#	modified by Kai Tang 
#	Jan 11,2011
# --------------------------------------------------------------------


# --------------------------------------------------------------------
#  ShoreMap extension to SHORE:
#  Identification of causal mutations using IlluminaGA2 sequencing data
# 
#  Written by Stephan Ossowski and Korbinian Schneeberger
#  --------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $referror = "";
my $marker;
my $consensus;
my $chrsizes;
my $window_step;
my $window_size_string;

my %CMD;
GetCom();



### Read in Ref Errors ------------------------------------------------------------------------
my %REFERROR = ();

if ($referror ne "") {
	open REFERRORFILE, $referror;
	while (my $line = <REFERRORFILE>) {
		my @a = split " ", $line;
		$REFERROR{$a[0]."#".$a[1]} = 1;
	}
	close REFERRORFILE;
}



### Read in Marker ----------------------------------------------------------------------------
my %MARKER = ();

open MARKERFILE, $marker or die "Cannot open marker file\n";
while (my $line=<MARKERFILE>) {
	my @a = split " ", $line;
	if (not defined($REFERROR{$a[0]."#".$a[1]})) {
		$MARKER{$a[0]."#".$a[1]} = $a[3];
	}
}
close MARKERFILE or die "Marker file won't close\n";

print STDERR "Read in marker file\n";



### Get the counts at the marker positions ----------------------------------------------------
my @SNPCHR = ();
my @SNPPOS = ();
my @REFFREQ = ();
my @ALLELEFREQ = ();

my $curr_chr = -1;


###########
#
# main change is here
#
###########


open CONS, $consensus or die "Cannot open consensus file\n";
while(my $line = <CONS>) {
	my @a = split " ", $line;
	my $chromosome = $a[0];
	my $position = $a[1];	
	my $coverage = $a[3];
	
	#	my $refbase = uc($a[44]);
	my $refbase = uc($a[2]);
	
	
	if ($chromosome ne $curr_chr) {
		print STDERR "Reading in chromosome: ", $chromosome, "\n";
		$curr_chr = $chromosome;
	}

	if ( defined($MARKER{$chromosome."#".$position} && 
	    ($refbase eq "A" || $refbase eq "C" || $refbase eq "G" || $refbase eq "T"))) {

		my $snpbase = uc($MARKER{$chromosome."#".$position});

		my $snp_count = 0;
		my $ref_count = 0;

		if ($coverage != 0) {

			if ($snpbase eq "A") {
				$snp_count = $a[4];
			}
			if ($snpbase eq "C") {
                                $snp_count = $a[5];
                        }
			if ($snpbase eq "G") {
                                $snp_count = $a[6];
                        }
			if ($snpbase eq "T") {
                                $snp_count = $a[7];
                        }

			if ($refbase eq "A") {
				$ref_count = $a[4];
			}
			if ($refbase eq "C") {
                                $ref_count = $a[5];
                        }
			if ($refbase eq "G") {
                                $ref_count = $a[6];
                        }
			if ($refbase eq "T") {
                                $ref_count = $a[7];
                        }
		}

	        push @SNPCHR, $chromosome;
	        push @SNPPOS, $position;
        	push @REFFREQ, $ref_count;
        	push @ALLELEFREQ, $snp_count;
	
		print $chromosome, "\t", $position, "\n";
	}
}

print STDERR "Read in consensus file\n";



### Create sliding windows --------------------------------------------------------------------

# Set windowsizes
my @window_sizes = split ",", $window_size_string;

for (my $w = 0; $w < @window_sizes; $w++) {
	my $window_size = $window_sizes[$w];

	print STDERR "Analysing window size: ", $window_size, "\n";

	my $output_string = "";

	my $chromosome = $SNPCHR[0];
	my $last_pos = $SNPPOS[0] - 1;
	my $marker_count = 0;
	my $ref_count = 0;
	my $allele_count = 0;

	my $j = 0;
	for (my $i = 0; $i < @SNPCHR; $i++) {

        	my $chr = $SNPCHR[$i];
	        my $pos = $SNPPOS[$i];
        	my $ref = $REFFREQ[$i];
	        my $mut = $ALLELEFREQ[$i];

		# delete old values
		
		
		#while ($SNPCHR[$j] != $chromosome) {
	     
		while ($SNPCHR[$j] ne $chromosome) {
					$marker_count = 0;
	                $ref_count = 0;
        	        $allele_count = 0;   
			$j++;
       		 }
	        $chromosome = $chr;
        	while($SNPPOS[$j] < $pos - $window_size) {
                	$marker_count--;
			if ($marker_count < 0) {
				$marker_count = 0;
			}
	                $ref_count -= $REFFREQ[$j];
			if ($ref_count < 0) {
				$ref_count = 0;
			}
        	        $allele_count -= $ALLELEFREQ[$j];
			if ($allele_count < 0) {
				$allele_count = 0;
			}
                	$j++;
        	}

		# add new value
		$marker_count++;
	        $ref_count += $REFFREQ[$i];
        	$allele_count += $ALLELEFREQ[$i];

		# Report if more than
		if ($pos > $last_pos + $window_step) {
                	if ($marker_count!= 0) {

				$output_string .= $chromosome."\t";
				$output_string .= int($pos-(($pos-$SNPPOS[$j])/2))."\t";
				$output_string .= $ref_count."\t".$allele_count."\t";

	                        if ($ref_count > $allele_count) {
        	                        if ($allele_count != 0) {
						$output_string .= $ref_count / $allele_count;
                        	        }
                                	else {
						$output_string .= $ref_count;	
	                                }
        	                }
                	        else {
                        	        if ($ref_count == $allele_count) {
						$output_string .= "0";
	                                }
        	                        else {
                	                        if ($ref_count != 0) {
							$output_string .= (-1) * ($allele_count / $ref_count);
                                	        }
                                        	else {
							$output_string .= "-".$allele_count;
	                                        }
        	                        }
                	        }
				$output_string .= "\n";
			}
                }
        }


	### Output the whole thing
	my $outputfile = "SHOREmap.winstep".$window_step.".winsize".$window_size.".txt";
	open OUT, "> ".$outputfile or die "Cannot open outputfile\n";
	print OUT $output_string;
	close OUT;

	### Call R for plotting
	my $pdffile = "SHOREmap.winstep".$window_step.".winsize".$window_size.".pdf";
	my $cmd = "R --slave --vanilla --args ".$chrsizes." ".$pdffile." ".$outputfile." ".$window_step." ".$window_size."  < /Users/tang58/scripts_all/perl_code/Kai_modified_SHOREmap/Kai_SHOREmap.R";
	print STDERR $cmd, "\n";
	system($cmd); 

}


### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

Mandatory:
--consensus    STRING   SHORE consensus file
--marker       STRING   Marker file
--chrsizes     STRING   Chromosome sizes file

Optional:
--referrors    STRING   Known reference sequence errors file
--windowstep   INT      Number of bp between datapoints. Default is 10000.
--windowsize   STRING   Comma separated windowsizes. 
\t\tDefault is: 50000,100000,150000,200000,250000,300000,350000,400000,450000,500000

See documentation for file formats.

SHOREmap is written by Korbinian Schneeberger and Stephan Ossowski.
Max Planck Institute for Developemental Biology, Tübingen, 2009.

\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "marker=s","consensus=s", "chrsizes=s", "windowsize=s", "windowstep=s", "referrors=s");

        die("Please specify marker file\n") unless defined($CMD{marker});
        die("Please specify consensus file\n") unless defined($CMD{consensus});
        die("Please specify chromosome sizes file\n") unless defined($CMD{chrsizes});

        $marker = $CMD{marker};
        $consensus = $CMD{consensus};
	$chrsizes = $CMD{chrsizes};

	if (defined($CMD{referrors})) {
		$referror = $CMD{referrors};
	}

	if (defined($CMD{windowstep})) {
		$window_step = $CMD{windowstep};
	}
	else {
		$window_step = 10000;
	}

	if (defined($CMD{windowsize})) {
                $window_size_string = $CMD{windowsize};
        }
        else {
                $window_size_string = "50000,100000,150000,200000,250000,300000,350000,400000,450000,500000";
        }
}


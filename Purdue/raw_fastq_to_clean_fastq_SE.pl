#!/usr/bin/perl -w
# Copyright (C)  2011-
# Program:			raw_fastq_to_clean_fastq_SE.pl
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		July 4, 2011
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	
# Description:		modified from Dr. Renyi Liu's (UCR) program
#**************************
# Version: 1.0 input fastq, output fastq
# 
#**************************
# e-mail:tangkai.ustc@gmail.com

# diff from orignal:
# orignal requries the length from 1 and 2 are the same

my $version="1.0";
use strict;

use Getopt::Long;
#use FindBin;
#use Bio::SeqIO;

my $debug = 0;
### Command line parameters -------------------------------------------------------------------
my $input	= "";	#input fastq file mate 1;
my $adaptor_P1	= "";	#mate 1 adaptor string;
my $minLen		= "";
my $phred = "";
my $output_pre = "";
my $outdir = "";
my %CMD;
my $poor_quality_char = "@";
GetCom();

my $ascii_cutoff = chr ($phred + ord($poor_quality_char));

print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";

die "wrong output dir" unless (-d $outdir);

my $P1_adaptor_1st_1nt = substr($adaptor_P1, 0, 1);
my $P1_adaptor_1st_2nt = substr($adaptor_P1, 0, 2);
my $P1_adaptor_1st_3nt = substr($adaptor_P1, 0, 3);
my $P1_adaptor_1st_4nt = substr($adaptor_P1, 0, 4);
my $P1_adaptor_1st_5nt = substr($adaptor_P1, 0, 5);
my $P1_adaptor_1st_6nt = substr($adaptor_P1, 0, 6);
my $P1_adaptor_1st_7nt = substr($adaptor_P1, 0, 7);
my $P1_adaptor_1st_8nt = substr($adaptor_P1, 0, 8);
my $P1_adaptor_1st_9nt = substr($adaptor_P1, 0, 9);
#my $P1_adaptor_1st_10nt = substr($adaptor_P1, 0, 10);
#my $P1_adaptor_1st_11nt = substr($adaptor_P1, 0, 11);
#my $P1_adaptor_1st_12nt = substr($adaptor_P1, 0, 12);
#my $P1_adaptor_1st_13nt = substr($adaptor_P1, 0, 13);

#my $P1_adaptor_1st_13nt_rev = reverse $P1_adaptor_1st_13nt;
my $P1_adaptor_1st_9nt_rev = reverse $P1_adaptor_1st_9nt;

my ($P1_trim_1, $P1_trim_2, $P1_trim_3, $P1_trim_4, $P1_trim_5, $P1_trim_6, $P1_trim_7,
	$P1_trim_8, $P1_trim_9_over, $P1_no_trim) = (0) x 10 ; # number of seqs subject to different trim
	
my $P1_no_short_seqs = 0;	

	

print STDERR "Handling files:\n$input\n";

#my $seqin	= Bio::SeqIO->new(-file=>$input, -format=>'fastq');
#my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fasta',width=>200);
############################# 1.1
#my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fastq',width=>200);
##############################

my $output = $output_pre."_clean.fastq"; 
 

print STDERR "output:\n", $output , "\n\n";
		
die "wrong input\n" unless (-e $input);

die "output exists\n" if (-e "$outdir/$output");
	
open (IN, $input) or die "cannot open $input:$!";

open (OUT, ">$outdir/$output") or die "cannot open $output: $!";

#for the trim low quality stat
my @P1_segment_hist;
my $P1_segment_sum   = 0;

#for the trim low quality stat mate2

my $single_P1 = 0;
my @S1_segment_hist;
my $S1_segment_sum = 0;


my $count = 0;

my $short = 0;

while(my $l1 = <IN>){
	$count++;
	if($count % 3000000 == 0){
		print STDERR "$count\tDONE\n";
	}
	my $P1_ID1 = $l1; 		
		
	chomp(my $P1_seq_string = <IN>);
	
	my $P1_ID2 = <IN>;
		
	chomp( my $P1_quality_string = <IN> );
		
		my $P1_original_length  = length $P1_seq_string;
		
		my $P1_cutoff_hit       =  0;
		my $P1_best_start_index =  0;
		my $P1_best_length      =  0;
		my $P1_current_start    =  0;
		#for P1
	#	for( my $i = 0; $i < $original_length; $i++ ){
		for( my $i = 0; $i < $P1_original_length; $i++ ){
				# if the quality score at this position is worse than the cutoff
				if( substr($P1_quality_string, $i, 1) le $ascii_cutoff ){
					$P1_cutoff_hit = 1;
					# determine length of good segment that just ended
					my $P1_current_segment_length = $i - $P1_current_start;
					# if this segment is the longest so far
					if( $P1_current_segment_length > $P1_best_length ){
						# store this segment as current best
						$P1_best_length      = $P1_current_segment_length;
						$P1_best_start_index = $P1_current_start;
					}
					# reset current start
					$P1_current_start = $i + 1;
			#	}elsif( $i == $original_length - 1){
				}elsif( $i == $P1_original_length - 1){
					# determine length of good segment that just ended
					my $P1_current_segment_length = ($i + 1) - $P1_current_start;
					# if this segment is the longest so far
					if( $P1_current_segment_length > $P1_best_length ){
						# store this segment as current best
						$P1_best_length = $P1_current_segment_length;
						$P1_best_start_index = $P1_current_start;
					}
				}
		}# for
			
		if( !$P1_cutoff_hit ){
			#$P1_best_length = $original_length;
			$P1_best_length = $P1_original_length;
		}
		$P1_segment_hist[ $P1_best_length ]++;
		$P1_segment_sum += $P1_best_length;
		
		if ($P1_best_length <= 0) {
      		$P1_seq_string = "N";
      		$P1_quality_string = $poor_quality_char;
    	} else {
     		$P1_seq_string = substr($P1_seq_string, $P1_best_start_index, $P1_best_length);
      		$P1_quality_string = substr($P1_quality_string, $P1_best_start_index, $P1_best_length);
    	}


##########################################################################################################
		my $P1_seqstr = $P1_seq_string;
		my $P1_new_seqstr = $P1_seqstr;
	
	if(length($P1_seqstr) < $minLen){
		#do nothing
	}
	else{	
		my $P1_flag_no	= 0;	my $P1_flag_ge9	= 0;
		if($P1_seqstr =~ /(\w*$P1_adaptor_1st_9nt)\w*/){			
			$P1_new_seqstr = $1;
			$P1_trim_9_over++;
			$P1_new_seqstr = reverse $P1_new_seqstr;
			if ($P1_new_seqstr =~ /\w*($P1_adaptor_1st_9nt_rev\w*)/){
				$P1_new_seqstr = $1;
			}
			$P1_new_seqstr = reverse $P1_new_seqstr;
			$P1_new_seqstr =~ s/$P1_adaptor_1st_9nt//;

		}elsif(substr($P1_seqstr, length($P1_seqstr) - 8) eq $P1_adaptor_1st_8nt){
			$P1_trim_8++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 8);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 7) eq $P1_adaptor_1st_7nt){
			$P1_trim_7++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 7);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 6) eq $P1_adaptor_1st_6nt){
			$P1_trim_6++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 6);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 5) eq $P1_adaptor_1st_5nt){
			$P1_trim_5++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 5);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 4) eq $P1_adaptor_1st_4nt){
			$P1_trim_4++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 4);
		
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 3) eq $P1_adaptor_1st_3nt){
			$P1_trim_3++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 3);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 2) eq $P1_adaptor_1st_2nt){
			$P1_trim_2++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 2);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 1) eq $P1_adaptor_1st_1nt){
			$P1_trim_1++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 1);
			
		}else{
			$P1_no_trim++;
			
		}

		if(length($P1_new_seqstr) < $minLen){
			$P1_no_short_seqs++;
			$P1_seq_string = $P1_new_seqstr; 
		}
		else{
			$P1_seq_string = $P1_new_seqstr ;
      		$P1_quality_string = substr($P1_quality_string, 0 , length($P1_new_seqstr));
		}	
	}
##########################################################################################################
	
##########################################################################################################
	if( length($P1_seq_string) >= $minLen ){
		
		my $temp_len_1 = length($P1_seq_string);
		
		
		print OUT $P1_ID1, $P1_seq_string, "\n+\n", $P1_quality_string, "\n";
		
		
		$single_P1++;
		
		$S1_segment_hist[$temp_len_1 ]++;
		$S1_segment_sum += $temp_len_1;
		
		
	}else{
		$short++;
	}
}		
	

summary();
sub_end_program();


############################################################################################################
######################                  summary
############################################################################################################
sub summary{
	
	my $mean_low_P1 = sprintf("%.1f",  $P1_segment_sum / $count);
	my $mean_ada_S1 = sprintf("%.1f",  $S1_segment_sum / $single_P1);
	
	my $median_low_P1 = find_median(\@P1_segment_hist, $count);
	my $median_ada_S1 = find_median(\@S1_segment_hist , $single_P1) ;

#	print "inputs:$input_P1\t$input_P2\ntotal input reads for each: $count\noutput:\n";
	print "inputs:$input\ntotal input reads for each: $count\noutput:\n";
	print "output pairs:$output\t$single_P1\n";

	print "discard too short pairs (<$minLen): $short\n\n";
	
	print "trim_low_qality for $input:\n";
	print "mean segment length = $mean_low_P1, median segment length = $median_low_P1\n" ;

	
	print "stat for $output:\n";
	print "mean segment length = $mean_ada_S1, median segment length = $median_ada_S1\n" ;

	
	

	
	print "for $input:\n";
	
	print join("\t",(">=9", "8", "7", "6", "5", "4", "3", "2", "1", "0",
	"short than $minLen bp")), "\n";
	print join("\t",($P1_trim_9_over, $P1_trim_8,
	$P1_trim_7, $P1_trim_6, $P1_trim_5, $P1_trim_4, $P1_trim_3, $P1_trim_2,
	$P1_trim_1, $P1_no_trim, $P1_no_short_seqs)),"\n";
		
}

############################################################################################################
######################                  find_median
############################################################################################################
sub find_median{
	my ($ref, $total) = @_;
	my $current_sum   = 0;
	my $current_index = 0;
	my $median_index1;
	my $median_index2;
	
	my $halfway_index = $total / 2;
	
	# while median_index1 and median_index2 are not defined
	while( !defined( $median_index1 ) || !defined( $median_index2 ) ){
		
		# add segment count to current sum for each segment length from array
		if( defined( $ref -> [ $current_index ] ) ){
			$current_sum += ( $ref -> [ $current_index ] ) ;
		}
		# if current sum of segment counts has surpassed halfway index
		if( $current_sum > $halfway_index ){
			# if median_index1 has not been defined, store current segment length
			if( !defined( $median_index1 ) ){
				$median_index1 = $current_index;
			}
			# if median_index2 has not been defined, store current segment length
			if( !defined( $median_index2 ) ){
				$median_index2 = $current_index;
			}
			# else if current sum of segment counts is exactly equal to the halfway index
		}elsif( $current_sum == $halfway_index	&& !defined( $median_index1 ) ){
			# store current segment length as median_index1
			$median_index1 = $current_index;
		}
		# loop through all possible segment lengths
		$current_index++;
	}
	# if number of segments is odd, store index2 as median segment length
	if( $total % 2 == 1){
		return $median_index1;
		# if number of segments is even, store average of index1 and index2 as median segment length
	}else{
		my $segment_median = sprintf( "%.0f", ( ( $median_index1 + $median_index2 ) / 2 ) );
		return $segment_median ;
	}
}
############################################################################################################
######################                  Read command line parameters
############################################################################################################
sub GetCom {
	#	my @usage = ("\nUsage: $0
	
	#	Mandatory:
	#--adaptor_P1	STRING	adapter_pair1
	#--adaptor_P2	STRING	adapter_pair2
	#--input_P1		STRING	fastq_file_pair1
	#--input_P2		STRING  fastq_file_pair2
	#--output_pre	STRING	prefix of output
	#--dir			STRING	output dir	
	#--minLen		Number	minium length of seq can be left after triming
	#--phred		    Phred score (between 0 and 40) at which base-calling error is considered too high
	#\n");
	
	my @usage = ("\nUsage: $0
	
	Mandatory:
	--adaptor	STRING	adapter
	--input		STRING	fastq_file_pair1
	--output_pre	STRING	prefix of output
	--dir			STRING	output dir	
	--minLen		Number	minium length of seq can be left after triming
	--phred		    Phred score (between 0 and 40) at which base-calling error is considered too high
	\n");
	
	
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "adaptor=s",   "input=s", "output_pre=s","minLen=n", "phred=n", "dir=s");
	
	# Mandatory params
	die("Please specify adaptor sequence\n") unless defined($CMD{adaptor});
	#	die("Please specify adaptor_P2 sequence\n") unless defined($CMD{adaptor_P2});
	
	die("Please specify input file\n") unless defined($CMD{input});
	#	die("Please specify input_P2 file\n") unless defined($CMD{input_P2});
	die("Please specify output dir\n") unless defined($CMD{dir});
	die("Please specify output prefix\n") unless defined($CMD{output_pre});
	die("Please specify minLen\n")	unless defined($CMD{minLen});
	die("Please specify phred score\n")	unless defined($CMD{phred});
	
	$adaptor_P1	= $CMD{adaptor};	#adaptor string;
	#	$adaptor_P2	= $CMD{adaptor_P2};	#adaptor string;
	$input	= $CMD{input};	#input fastq file;
	#	$input_P2	= $CMD{input_P2};	#input fastq file;
	$output_pre = $CMD{output_pre};	#output fasta file;
	$minLen		= $CMD{minLen};
	$phred 		= $CMD{phred};
	$outdir		= $CMD{dir};

	print STDERR "adaptor:$adaptor_P1\n", "input:$input、n", 
					"minium length is $minLen\nmin_phred:$phred\noutdir:$outdir\n";
}


############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);
	
}
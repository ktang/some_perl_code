#!/usr/bin/perl -w
# Copyright (C)  2011-
# Program:			raw_fastq_to_clean_fastq_v1.0.pl
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

my $version="1.0";
use strict;

use Getopt::Long;
use FindBin;
use Bio::SeqIO;

my $debug = 0;
### Command line parameters -------------------------------------------------------------------
my $input_P1	= "";	#input fastq file mate 1;
my $input_P2	= "";	#input fastq file mate 1;
my $adaptor_P1	= "";	#mate 1 adaptor string;
my $adaptor_P2	= "";	#mate 2 adaptor string;
my $minLen		= "";
my $phred = "";
my $output_pre = "";
my $outdir = "";
my %CMD;
my $poor_quality_char = "@";
GetCom();

my $ascii_cutoff = chr ($phred + 64);

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
my $P1_adaptor_1st_10nt = substr($adaptor_P1, 0, 10);
my $P1_adaptor_1st_11nt = substr($adaptor_P1, 0, 11);
my $P1_adaptor_1st_12nt = substr($adaptor_P1, 0, 12);
my $P1_adaptor_1st_13nt = substr($adaptor_P1, 0, 13);

my $P1_adaptor_1st_13nt_rev = reverse $P1_adaptor_1st_13nt;

my ($P1_trim_1,$P1_trim_2, $P1_trim_3, $P1_trim_4, $P1_trim_5, $P1_trim_6, $P1_trim_7,
	$P1_trim_8, $P1_trim_9, $P1_trim_10, $P1_trim_11, $P1_trim_12, $P1_trim_13_over, $P1_no_trim) =
	(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); # number of seqs subject to different trim
my $P1_no_short_seqs = 0;	

my $P2_adaptor_1st_1nt = substr($adaptor_P2, 0, 1);
my $P2_adaptor_1st_2nt = substr($adaptor_P2, 0, 2);
my $P2_adaptor_1st_3nt = substr($adaptor_P2, 0, 3);
my $P2_adaptor_1st_4nt = substr($adaptor_P2, 0, 4);
my $P2_adaptor_1st_5nt = substr($adaptor_P2, 0, 5);
my $P2_adaptor_1st_6nt = substr($adaptor_P2, 0, 6);
my $P2_adaptor_1st_7nt = substr($adaptor_P2, 0, 7);
my $P2_adaptor_1st_8nt = substr($adaptor_P2, 0, 8);
my $P2_adaptor_1st_9nt = substr($adaptor_P2, 0, 9);
my $P2_adaptor_1st_10nt = substr($adaptor_P2, 0, 10);
my $P2_adaptor_1st_11nt = substr($adaptor_P2, 0, 11);
my $P2_adaptor_1st_12nt = substr($adaptor_P2, 0, 12);
my $P2_adaptor_1st_13nt = substr($adaptor_P2, 0, 13);

my $P2_adaptor_1st_13nt_rev = reverse $P2_adaptor_1st_13nt;

my ($P2_trim_1,$P2_trim_2, $P2_trim_3, $P2_trim_4, $P2_trim_5, $P2_trim_6, $P2_trim_7,
	$P2_trim_8, $P2_trim_9, $P2_trim_10, $P2_trim_11, $P2_trim_12, $P2_trim_13_over, $P2_no_trim) =
	(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); # number of seqs subject to different trim
my $P2_no_short_seqs = 0;	

print STDERR "Handling files:\n$input_P1\n$input_P2\n";

#my $seqin	= Bio::SeqIO->new(-file=>$input, -format=>'fastq');
#my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fasta',width=>200);
############################# 1.1
#my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fastq',width=>200);
##############################

my $output_P1 = $output_pre."_1_pairs.fastq"; 
my $output_P2 = $output_pre."_2_pairs.fastq"; 
my $output_S1 = $output_pre."_1_single.fastq"; 
my $output_S2 = $output_pre."_2_single.fastq"; 

print STDERR "output files:\t", join("\n", ($output_P1, 
		$output_P2, $output_S1, $output_S2)), "\n";
		
die "wrong input\n" unless (-e $input_P1 and -e $input_P2);

die "output exists\n" unless (!(-e "$outdir/$output_P1") and !(-e "$outdir/$output_P2") 
	and !(-e "$outdir/$output_S1") and !(-e "$outdir/$output_S2"));
	
open (IN1, $input_P1) or die "cannot open $input_P1:$!";
open (IN2, $input_P2) or die "cannot open $input_P2):$!";
open (PAIR1, ">$outdir/$output_P1") or die "cannot open $output_P1:$!";
open (PAIR2, ">$outdir/$output_P2") or die "cannot open $output_P2:$!";
open (SI1, ">$outdir/$output_S1") or die "cannot open $output_S1:$!";
open (SI2, ">$outdir/$output_S2") or die "cannot open $output_S2:$!";

#for the trim low quality stat
my @P1_segment_hist;
my $P1_segment_sum   = 0;

#for the trim low quality stat mate2
my @P2_segment_hist;
my $P2_segment_sum   = 0;

my $paired_no = 0;
my @P1_paired_segment_hist;
my $P1_paired_segment_sum = 0;

my @P2_paired_segment_hist;
my $P2_paired_segment_sum = 0;

my $single_P1 = 0;
my @S1_segment_hist;
my $S1_segment_sum = 0;

my $single_P2 = 0;
my @S2_segment_hist;
my $S2_segment_sum = 0;

my $count = 0;

my $both_short = 0;

while(my $l1 = <IN1>, my $l2 = <IN2>){
	$count++;
	if($count % 3000000 == 0){
		print STDERR "$count\tDONE\n";
	}
		my $P1_ID1 = $l1; 		
		my $P2_ID1 = $l2;
		
		chomp(my $P1_seq_string = <IN1>);
		chomp(my $P2_seq_string = <IN2>);
		
		my $P1_ID2 = <IN1>;
		my $P2_ID2 = <IN2>;
		
		chomp( my $P1_quality_string = <IN1> );
		chomp( my $P2_quality_string = <IN2> );
		
		my $P1_original_length  = length $P1_seq_string;
		my $P2_original_length  = length $P2_seq_string;
		if ($P1_original_length != $P2_original_length){
			print STDERR "wrong length\n$P1_ID1\n$P1_seq_string\n$P2_ID1\n$P2_seq_string\n";
			die;
		}
		my $original_length = $P1_original_length;
###################################################################################
		my $P1_cutoff_hit       =  0;
		my $P1_best_start_index =  0;
		my $P1_best_length      =  0;
		my $P1_current_start    =  0;
		#for P1
		for( my $i = 0; $i < $original_length; $i++ ){
		
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
				}elsif( $i == $original_length - 1){
					# determine length of good segment that just ended
					my $P1_current_segment_length = ($i + 1) - $P1_current_start;
					# if this segment is the longest so far
					if( $P1_current_segment_length > $P1_best_length ){
						# store this segment as current best
						$P1_best_length = $P1_current_segment_length;
						$P1_best_start_index = $P1_current_start;
					}
				}
		}
			
		if( !$P1_cutoff_hit ){
			$P1_best_length = $original_length;
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
###################################################################################
		my $P2_cutoff_hit       =  0;
		my $P2_best_start_index =  0;
		my $P2_best_length      =  0;
		my $P2_current_start    =  0;
		for( my $i = 0; $i < $original_length; $i++ ){
		
				if( substr($P2_quality_string, $i, 1) le $ascii_cutoff ){
					$P2_cutoff_hit = 1;
					my $P2_current_segment_length = $i - $P2_current_start;
					if( $P2_current_segment_length > $P2_best_length ){
						$P2_best_length      = $P2_current_segment_length;
						$P2_best_start_index = $P2_current_start;
					}
					$P2_current_start = $i + 1;
				}elsif( $i == $original_length - 1){
					my $P2_current_segment_length = ($i + 1) - $P2_current_start;
					if( $P2_current_segment_length > $P2_best_length ){
						$P2_best_length = $P2_current_segment_length;
						$P2_best_start_index = $P2_current_start;
					}
				}
		}

		if( !$P2_cutoff_hit ){
			$P2_best_length = $original_length;
		}
	
		$P2_segment_hist[ $P2_best_length ]++;
		$P2_segment_sum += $P2_best_length;
		
		if ($P2_best_length <= 0) {
      		$P2_seq_string = "N";
      		$P2_quality_string = $poor_quality_char;
    	} else {
     		$P2_seq_string = substr($P2_seq_string, $P2_best_start_index, $P2_best_length);
      		$P2_quality_string = substr($P2_quality_string, $P2_best_start_index, $P2_best_length);
    	}
##########################################################################################################
		my $P1_seqstr = $P1_seq_string;
		my $P1_new_seqstr = $P1_seqstr;
	if(length($P1_seqstr) < $minLen){
		#do nothing
	}
	else{	
		my $P1_flag_no	= 0;	my $P1_flag_ge13	= 0;
		if($P1_seqstr =~ /(\w*$P1_adaptor_1st_13nt)\w*/){			
			$P1_new_seqstr = $1;
			$P1_trim_13_over++;
			$P1_new_seqstr = reverse $P1_new_seqstr;
			if ($P1_new_seqstr =~ /\w*($P1_adaptor_1st_13nt_rev\w*)/){
				$P1_new_seqstr = $1;
			}
			$P1_new_seqstr = reverse $P1_new_seqstr;
			$P1_new_seqstr =~ s/$P1_adaptor_1st_13nt//;

		}elsif(substr($P1_seqstr, length($P1_seqstr) - 12) eq $P1_adaptor_1st_12nt){
			$P1_trim_12++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 12);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 11) eq $P1_adaptor_1st_11nt){
			$P1_trim_11++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 11);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 10) eq $P1_adaptor_1st_10nt){
			$P1_trim_10++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 10);
			
		}elsif(substr($P1_seqstr, length($P1_seqstr) - 9) eq $P1_adaptor_1st_9nt){
			$P1_trim_9++;
			$P1_new_seqstr = substr($P1_seqstr, 0, length($P1_seqstr) - 9);
			
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
	my $P2_seqstr = $P2_seq_string;
	my $P2_new_seqstr = $P2_seqstr;
	if(length($P2_seqstr) < $minLen){
		#do nothing
	}	else{	
		if($P2_seqstr =~ /(\w*$P2_adaptor_1st_13nt)\w*/){			
			$P2_new_seqstr = $1;

			$P2_trim_13_over++;

			$P2_new_seqstr = reverse $P2_new_seqstr;
			if ($P2_new_seqstr =~ /\w*($P2_adaptor_1st_13nt_rev\w*)/){
				$P2_new_seqstr = $1;
			}
			$P2_new_seqstr = reverse $P2_new_seqstr;
			$P2_new_seqstr =~ s/$P2_adaptor_1st_13nt//;
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 12) eq $P2_adaptor_1st_12nt){
			$P2_trim_12++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 12);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 11) eq $P2_adaptor_1st_11nt){
			$P2_trim_11++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 11);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 10) eq $P2_adaptor_1st_10nt){
			$P2_trim_10++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 10);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 9) eq $P2_adaptor_1st_9nt){
			$P2_trim_9++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 9);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 8) eq $P2_adaptor_1st_8nt){
			$P2_trim_8++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 8);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 7) eq $P2_adaptor_1st_7nt){
			$P2_trim_7++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 7);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 6) eq $P2_adaptor_1st_6nt){
			$P2_trim_6++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 6);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 5) eq $P2_adaptor_1st_5nt){
			$P2_trim_5++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 5);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 4) eq $P2_adaptor_1st_4nt){
			$P2_trim_4++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 4);
			
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 3) eq $P2_adaptor_1st_3nt){
			$P2_trim_3++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 3);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 2) eq $P2_adaptor_1st_2nt){
			$P2_trim_2++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 2);
			
		}elsif(substr($P2_seqstr, length($P2_seqstr) - 1) eq $P2_adaptor_1st_1nt){
			$P2_trim_1++;
			$P2_new_seqstr = substr($P2_seqstr, 0, length($P2_seqstr) - 1);
			
		}else{
			$P2_no_trim++;	
			
		}
		
		if(length($P2_new_seqstr) < $minLen){
			$P2_no_short_seqs++;
			$P2_seq_string = $P2_new_seqstr;

		}
		else{
			$P2_seq_string = $P2_new_seqstr;
      		$P2_quality_string = substr($P2_quality_string, 0 , length($P2_new_seqstr));
		}	
	}
##########################################################################################################
	if( length($P1_seq_string) >= $minLen and length($P2_seq_string) >= $minLen ){
		
		my $temp_len_1 = length($P1_seq_string);
		my $temp_len_2 = length($P2_seq_string);
		
		print PAIR1 $P1_ID1, $P1_seq_string, "\n+\n", $P1_quality_string, "\n";
		print PAIR2 $P2_ID1, $P2_seq_string, "\n+\n", $P2_quality_string, "\n";
		
		$paired_no++;
		
		$P1_paired_segment_hist[$temp_len_1 ]++;
		$P1_paired_segment_sum += $temp_len_1;
		
		$P2_paired_segment_hist[ $temp_len_2 ]++;
		$P2_paired_segment_sum += $temp_len_2;
	}elsif(length($P1_seq_string) >= $minLen){
		my $temp = length($P1_seq_string);
		
		print SI1 $P1_ID1, $P1_seq_string, "\n+\n", $P1_quality_string, "\n";
		$single_P1++;
		
		$S1_segment_hist[$temp]++;
		$S1_segment_sum += $temp;
	}elsif(length($P2_seq_string) >= $minLen){
		my $temp = length($P2_seq_string);

		print SI2 $P2_ID1, $P2_seq_string, "\n+\n", $P2_quality_string, "\n";
		
		$single_P2++;
		
		$S2_segment_hist[$temp]++;
		$S2_segment_sum += $temp;
	}else{
		$both_short++;
	}
}		
	

summary();
sub_end_program();


############################################################################################################
######################                  summary
############################################################################################################
sub summary{
	
	my $mean_low_P1 = sprintf("%.1f",  $P1_segment_sum / $count);
	my $mean_low_P2 = sprintf("%.1f",  $P2_segment_sum / $count);

	my $mean_ada_P1 = sprintf("%.1f",  $P1_paired_segment_sum / $paired_no);
	my $mean_ada_P2 = sprintf("%.1f",  $P2_paired_segment_sum / $paired_no);
	my $mean_ada_S1 = sprintf("%.1f",  $S1_segment_sum / $single_P1);
	my $mean_ada_S2 = sprintf("%.1f",  $S2_segment_sum / $single_P2);
	
	my $median_low_P1 = find_median(\@P1_segment_hist, $count) ;
	my $median_low_P2 = find_median(\@P2_segment_hist, $count) ;
	my $median_ada_P1 = find_median(\@P1_paired_segment_hist , $paired_no) ;
	my $median_ada_P2 = find_median(\@P2_paired_segment_hist , $paired_no) ;
	my $median_ada_S1 = find_median(\@S1_segment_hist , $single_P1) ;
	my $median_ada_S2 = find_median(\@S2_segment_hist , $single_P2) ;

	print "inputs:$input_P1\t$input_P2\ntotal input reads for each: $count\noutput:\n";
	print "output pairs:$output_P1\t$output_P2\npaired:$paired_no\n";
	print "single in mate1($output_S1): $single_P1\n";
	print "single in mate2($output_S2): $single_P2\n";
	print "discard too short pairs (<$minLen): $both_short\n\n";
	
	print "trim_low_qality for $input_P1:\n";
	print "mean segment length = $mean_low_P1, median segment length = $median_low_P1\n" ;

	print "trim_low_qality for $input_P2:\n";
	print "mean segment length = $mean_low_P2, median segment length = $median_low_P2\n\n" ;
	
	print "stat for $output_P1:\n";
	print "mean segment length = $mean_ada_P1, median segment length = $median_ada_P1\n" ;

	print "stat for $output_P2:\n";
	print "mean segment length = $mean_ada_P2, median segment length = $median_ada_P2\n" ;
	
	print "stat for $output_S1:\n";
	print "mean segment length = $mean_ada_S1, median segment length = $median_ada_S1\n" ;
	
	print "stat for $output_S2:\n";
	print "mean segment length = $mean_ada_S2, median segment length = $median_ada_S2\n\n" ;
	

	
	print "for $input_P1:\n";
	
	print join("\t",(">=13","12","11","10","9", "8", "7", "6", "5", "4", "3", "2", "1", "0",
	"short than $minLen bp")), "\n";
	print join("\t",($P1_trim_13_over,$P1_trim_12,$P1_trim_11,$P1_trim_10,$P1_trim_9, $P1_trim_8,
	$P1_trim_7, $P1_trim_6, $P1_trim_5, $P1_trim_4, $P1_trim_3, $P1_trim_2,
	$P1_trim_1, $P1_no_trim, $P1_no_short_seqs)),"\n";
	
	print "for $input_P2:\n";
	
	print join("\t",(">=13","12","11","10","9", "8", "7", "6", "5", "4", "3", "2", "1", "0",
	"short than $minLen bp")), "\n";
	print join("\t",($P2_trim_13_over,$P2_trim_12,$P2_trim_11,$P2_trim_10,$P2_trim_9, $P2_trim_8,
	$P2_trim_7, $P2_trim_6, $P2_trim_5, $P2_trim_4, $P2_trim_3, $P2_trim_2,
	$P2_trim_1, $P2_no_trim, $P2_no_short_seqs)),"\n";
	
	
	
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
  my @usage = ("\nUsage: $0

Mandatory:
--adaptor_P1	STRING	adapter_pair1
--adaptor_P2	STRING	adapter_pair2
--input_P1		STRING	fastq_file_pair1
--input_P2		STRING  fastq_file_pair2
--output_pre	STRING	prefix of output
--dir			STRING	output dir	
--minLen		Number	minium length of seq can be left after triming
--phred		    Phred score (between 0 and 40) at which base-calling error is considered too high
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "adaptor_P1=s",  "adaptor_P2=s", "input_P1=s", "input_P2=s",
	 "output_pre=s","minLen=n", "phred=n", "dir=s");
	
	# Mandatory params
	die("Please specify adaptor_P1 sequence\n") unless defined($CMD{adaptor_P1});
	die("Please specify adaptor_P2 sequence\n") unless defined($CMD{adaptor_P2});
	
	die("Please specify input_P1 file\n") unless defined($CMD{input_P1});
	die("Please specify input_P2 file\n") unless defined($CMD{input_P2});
	die("Please specify output dir\n") unless defined($CMD{dir});
	die("Please specify output prefix\n") unless defined($CMD{output_pre});
	die("Please specify minLen\n")	unless defined($CMD{minLen});
	die("Please specify phred score\n")	unless defined($CMD{phred});
	
	$adaptor_P1	= $CMD{adaptor_P1};	#adaptor string;
	$adaptor_P2	= $CMD{adaptor_P2};	#adaptor string;
	$input_P1	= $CMD{input_P1};	#input fastq file;
	$input_P2	= $CMD{input_P2};	#input fastq file;
	$output_pre = $CMD{output_pre};	#output fasta file;
	$minLen		= $CMD{minLen};
	$phred 		= $CMD{phred};
	$outdir		= $CMD{dir};

	print STDERR "adaptor_P1:$adaptor_P1\nadaptor_P2:$adaptor_P2\n", "input_P1:$input_P1\ninput_P2:$input_P2\n", 
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
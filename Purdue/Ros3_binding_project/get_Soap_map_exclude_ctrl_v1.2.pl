#!/usr/bin/perl -w
# Copyright (C)  	2011-
# Program:			get_Soap_map_exclude_ctrl.pl
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		Mar 21, 2011
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	
# Description:		input two soap outcome,one is control and the other is not 
#					get specific position which is in experiment not in control
#					
#**************************
# Version: 1.0	input two soap outcome,one is control and the other is not 
#				get specific position which is in experiment not in control
# Version: 1.1  not use a window, use exact postition(every 10 base)
# Version: 1.2  not use a window, use exact postition(record every position)
#**************************
# e-mail:tangkai.ustc@gmail.com

my $version="1.2";
use strict;
use warnings;

use Getopt::Long;


my $debug = 0;
### Command line parameters -------------------------------------------------------------------
my $input_ctrl	= "";	#input soap outcome of control
my $input_pos = ""; #input soap outcome of positive experiment.
my $output	= "";	#output fasta file;

#my $win_size = 100;
my %CMD;


print STDERR ("\n==================| $0 start |==========================================\n");
my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";

GetCom();

if (-e $output) {die "$output exists!!!\n"}
open(CTRL, $input_ctrl) or die "cannot open $input_ctrl:$!";

open(IN, $input_pos) or die "cannot open $input_pos:$!";

open (OUT, ">$output") or die "cannot open $output:$!";


my %ctrls;

while(<CTRL>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[7];
	my $pos = $a[8];
	my $end = $a[8] + $a[5] - 1;
#	my $win = int( ($pos - 1) / $win_size ) + 1;
	for my $i ($pos..$end){
		$ctrls{$chr} -> [$i] = 1;
	}
}

while(<IN>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[7];
	my $pos = $a[8];
	my $end = $a[8] + $a[5] - 1;
#	my $win = int( ($pos - 1) / $win_size ) + 1;
#	if ( !(defined $ctrls{$chr} -> [$win]) ){
#		print OUT $_,"\n";
#	my $itr = int($a[5] / 10);
	my $flag = 0;
#	for my $i (0..$itr){
#		my $index = $pos + $i * 10;
#		if (defined $ctrls{$chr}->[$index]){
#			$flag = 1;
#			last;
#		}
#	}

	for my $i ($pos..$end){
		if (defined $ctrls{$chr}->[$i]){
			$flag = 1;
			last;	
		}
	}	
	
	if ($flag == 0 and !(defined $ctrls{$chr}->[$end])){
		print OUT $_,"\n";
	}
	
#	}
}

sub_end_program();
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

############################################################################################################
######################                  Read command line parameters
############################################################################################################
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory:
--ctrl_input		STRING
--treatment_input		STRING
--output		STRING
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "ctrl_input=s", "treatment_input=s", "output=s");
	
	# Mandatory params
	die("Please specify ctrl_input file\n") unless defined($CMD{ctrl_input});
	die("Please specify treatment_input file\n") unless defined($CMD{treatment_input});
	die("Please specify output file\n") unless defined($CMD{output});
	$input_ctrl	= $CMD{ctrl_input};	#input fastq file;
	$output= $CMD{output};	#output fasta file;
	$input_pos = $CMD{treatment_input}
}


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

my $usage = "$0 <ctrl> <ip1> <out1> <ip2> <out2>";

die $usage unless (@ARGV == 5);

my ($ctrl, $ip1, $out1 , $ip2, $out2) = @ARGV[0..4];

print STDERR ("\n==================| $0 start |==========================================\n");
my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";


if (-e $out1) {die "$out1 exists!!!\n"}
if (-e $out2) {die "$out2 exists!!!\n"}
 
open(CTRL, $ctrl) or die "cannot open $ctrl:$!";

open(IP1, $ip1) or die "cannot open $ip1:$!";
open(IP2, $ip2) or die "cannot open $ip2:$!";

open (OUT1, ">$out1") or die "cannot open $out1:$!";
open (OUT2, ">$out2") or die "cannot open $out2:$!";

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



while(<IP1>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[7];
	my $pos = $a[8];
	my $end = $a[8] + $a[5] - 1;
	my $flag = 0;

	for my $i ($pos..$end){
		if (defined $ctrls{$chr}->[$i]){
			$flag = 1;
			last;	
		}
	}	
	
	if ($flag == 0 and !(defined $ctrls{$chr}->[$end])){
		print OUT1 $_,"\n";
	}
}

while(<IP2>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[7];
	my $pos = $a[8];
	my $end = $a[8] + $a[5] - 1;
	my $flag = 0;

	for my $i ($pos..$end){
		if (defined $ctrls{$chr}->[$i]){
			$flag = 1;
			last;	
		}
	}	
	
	if ($flag == 0 and !(defined $ctrls{$chr}->[$end])){
		print OUT2 $_,"\n";
	}
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


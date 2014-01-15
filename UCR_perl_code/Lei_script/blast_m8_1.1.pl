#!/usr/bin/perl
# Copyright (c)  2008-
# Program:			blast_m8
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2008.08.28
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2008.08.28
# Description:	analysis blast result with parameter m=8.
#**************************
# Version: 1.1	fix the output
#**************************
# e-mail:highlei@gmail.com

my $version="1.1";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:l:b:c:o:r:t:");
my $infile		= $opt_i;
my $lenfile		= (defined $opt_l) ? $opt_l : "";
my $batchfile	= $opt_b;
my $batchparam	= $opt_c;
my $output		= (defined $opt_o) ? $opt_o : $start;
my $lenRatio	= (defined $opt_r) ? $opt_r : 0.7;
my $identity	= (defined $opt_t) ? $opt_t : 95;

if ($opt_h){
	usage();
}

sub numerically {$a <=> $b}

use FileHandle;
use strict;



my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo);
my (%query,%len);
my $key="";


#===========================================================================================================
#====================                  main
#===========================================================================================================
my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %35s : %-25s","input blast result (m=8) file",$infile);#%45s
	printf STDERR ("\n l  %35s : %-25s","input length file",$lenfile);
	printf STDERR ("\n o  %35s : %-25s","output files prefix name",$output);
	printf STDERR ("\n r  %35s : %-25s","len threshold=align_len/query_len",$lenRatio);
	printf STDERR ("\n t  %35s : %-25s","identity threshold\%",$identity);
	printf STDERR ("\n x  %35s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n------------------------------------------------------------\n\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input the blast (m=8) result file:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "l") {print STDERR "please input the length file:\n"; $lenfile	= <STDIN>;$lenfile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "o") {print STDERR "please input output files prefix name:\n";$output	= <STDIN>;$output	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "r") {print STDERR "please input length threshold=align_len/query_len:\n";$lenRatio	= <STDIN>;$lenRatio	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "t") {print STDERR "please input identity threshold\%:\n";$identity	= <STDIN>;$identity	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(1);}
}

#$infile	=~s/&+/ /g;
#$batchfile	=~s/\"$//g;	$batchfile	=~s/^\"//g;
#$batchparam	=~s/\"$//g;	$batchparam	=~s/^\"//g;
#$output	=~s/\"$//g;		$output	=~s/^\"//g;
#$opparam	=~s/\"$//g;	$opparam	=~s/^\"//g;
#$param	=~s/\"$//g;		$param	=~s/^\"//g;
#
#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";

############################################ read file ######################################################

if ($lenfile ne "") {
	$file = new FileHandle ("$lenfile") || die("Cannot open file: $lenfile");

	while(<$file>)
	{
		if ($_ =~/^(\S+)\s+(\d+)/) {
			if (exists($len{$1})) {
				die("there are something wrong in $lenfile:$1\n");
			} else {
				$len{$1}	= $2;
			}
		} else {

		}
	}

	close $file || die;
	print	STDERR "read length file $lenfile completely!\n";
}

$k1	= 0;
$file = new FileHandle ("$infile") || die("Cannot open file: $infile");

while(<$file>)
{
	if ($_ =~/^(\S+)\s+(\S+)\s+([\d|\.]+)\s+(\d+)/) {
		$i	= $_;#	$i	=~s/[\s|\r|\n|\t]+$//g;
#		if (exists($query{$1}) && $query{$1}->[0]==1) {

#		} else {
			if ($lenfile ne "") {
				if ($4/$len{$1} < $lenRatio) {
					if (!exists($query{$1})) {
		#				$query{$1}->[0]	= 0;
						$k1++;
					}
					next;
				}
			}
			if ($3 < $identity) {
				if (!exists($query{$1})) {
		#			$query{$1}->[0]	= 0;
					$k1++
				}
				next;
			}
			if (!exists($query{$1})) {
				$k1++
			}
#			$query{$1}{"num"} = 1;
			if (exists($query{$1}{"num"}) ) {
				$query{$1}{"num"}++;
				$query{$1}{$query{$1}{"num"}}	= $i;
			} else {
				$query{$1}{"num"}	= 1;
				$query{$1}{$query{$1}{"num"}}	= $i;
			}
#			$query{$1}->[1]	= $i;
#		}
	} else {

	}
}

close $file || die;

print	STDERR "read file $infile completely,element num=$k1!\n";

############################################ output ######################################################
$i	= $output . "exist.lis";
open(EXIST, ">$i") || die("Can not open file: $i\n");
$j	= $output . "non_exist.lis";
#open(NONEXIST, ">$j") || die("Can not open file: $j\n");

$k1	= 0;	$k2	= 0;
foreach $key (sort keys(%query)) {
	foreach $k (sort numerically keys %{$query{$key}}) {
		if ($k	eq "num") {
			next;
		}
		print EXIST $query{$key}{$k};#,"\n";
		$k1++;
	}
#	if ($query{$key}->[0] == 1) {
#		print EXIST	$query{$key}->[1],"\n";
#		$k1++;
#	} else {
##		print NONEXIST	$key,"\n";
#		$k2++;
#	}
}
close (EXIST);

#close (NONEXIST);

print STDERR "exist: $k1\n";


sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################


sub usage
{
	print "Contact : Gaolei <highlei\@gmail.com>";
	print "\nProgram : $0\nVersion: $version\n";
	print "Usage:\n	$0 \n";
	print "-i	input the blast (m=8) result file\n";
	print "			eg: blast_m8.blastn\n";
	print "-l	input the length file\n";
	print "			eg: wgs_1.len\n";
#	print "-c	input batchfiles parameter\n";
#	print "			eg: input\n";
	print "-o	input output files prefix name\n";
	print "			eg: wgs_1_\n";
#	print "-m	output two gene files\n";
#	print "			eg: -m 1;default: donot output\n";
	print "-r	length threshold=align_len/query_len\n";
	print "			eg: -r 0.6; default: 0.5\n";
	print "-t	input identity threshold\%\n";
	print "			eg: -t 95; default: 90\n";
	print "-h	display this lines\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
	print "$0 -i blast_m8.blastn -t 90\n";
	print "$0 -i blast_m8.blastn -l wgs_1.len -r 0.5\n";
	print"\n============================================================\n\n";

    exit(1);
}
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #时间子程序
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n============================================================\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("............................................................\n");
#	exit(1);
	return;
}

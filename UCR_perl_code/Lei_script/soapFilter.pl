#!/usr/bin/perl
# Copyright (c)  2010-
# Program:			soap_filter
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2010.12.08
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2010.12.08
# Description:	read the soap result to filter
#**************************
# Version: 1.0
#**************************
# e-mail:highlei@gmail.com

my $version="1.0";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:p:c:f:b:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $posfile		= (defined $opt_p) ? $opt_p : "";
my $col			= (defined $opt_c) ? $opt_c : "8,9";
my $fore		= (defined $opt_f) ? $opt_f : 500;
my $back		= (defined $opt_b) ? $opt_b : 500;


if ($opt_h || $infile eq ""){
	usage();
}


sub numerically{$a<=>$b};

use FileHandle;
use strict;



my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$k4,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@genome,@gnmName,@gnmLen);
my (%gnm,%seg,%strand);
my $key="";
my ($endLen,$minLen,$addLen,$maxLen,$foldFileNum,$soap_M);
my ($foldFile,$foldFileOut,$soapOut);
my ($seq,$pos);
#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %55s : %-25s","input soap result",$infile);#%45s
	printf STDERR ("\n p  %55s : %-25s","input position file",$posfile);
	printf STDERR ("\n c  %55s : %-25s","input column",$col);
	printf STDERR ("\n f  %55s : %-25s","input fore length",$fore);
	printf STDERR ("\n b  %55s : %-25s","input back length",$back);#%45s
#	if($zero==1) {printf STDERR ("\n z  %45s : %-25s","output coverage region?","1");}
#	elsif ($zero == 0) {printf STDERR ("\n z  %45s : %-25s","concise output","0");}
#	else {printf STDERR ("\n z  %45s : %-25s","output all genome","2");}
	printf STDERR ("\n x  %55s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n------------------------------------------------------------\n\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input the soap result:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "p") {print STDERR "please input the position file:\n"; $posfile	= <STDIN>;$posfile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "c") {print STDERR "please input column:\n"; $col	= <STDIN>;$col	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input fore length\n";$fore	= <STDIN>;$fore	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "b") {print STDERR "please input back length\n";$fore	= <STDIN>;$fore	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(0);}
}

if ($col =~ /(\d+)\,(\d+)/) {
	$seq	= $1-1;	$pos	= $2-1;
	print STDERR "seq=$seq\tpos=$pos\n\n";
}

############################################ read file ######################################################

$k1	= 0;	$k3	= 0;
$file = new FileHandle ("$posfile") || die("Cannot open file: $posfile");

while(<$file>)
{
	@bufi	= split(/\t+/,$_);
	for ($k = $bufi[$pos] - $fore; $k < $bufi[$pos]+$back ;$k++) {
		$gnm{$bufi[$seq]}{$k}	= 1;	$k3++;
	}
#	$gnm{$bufi[$seq]}{$bufi[$pos]}	= 1;
	$k1++;
	splice(@bufi,0);
}

close $file || die;

print	STDERR "\nLoad pos $posfile OK\t$k1,$k3\n";
############################################ read soap files ######################################################
$k1	= 0;	$k3	= 0;
$file = new FileHandle ("$infile") || die("Cannot open file: $infile");

while(<$file>)
{
	if ($_=~/^\S+\s+\S+\s+\S+\s+\d+\s+\S+\s+\d+\s+[+|-]+\s+(\S+)\s+(\d+)\s+\d+\s+/) {
		$i	= $1;	$j	= $2;	$k3++;
		if (exists($gnm{$i}{$j})) {
			print $_;	$k1++;
			next;
		}
#		for ($k2 = $j - $fore; $k2 < $j+$back ;$k2++) {
#			if (exists($gnm{$i}{$k2})) {
#				print $_;	$k1++;
#				last;
#			}
#		}
	}
}

close $file || die;

print	STDERR "\nLoad soap file $infile OK\t$k1,$k3\n";


############################################ output ######################################################

sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################


sub usage
{
	print "Program :\t$0\n";
	print "Version :\t$version\n";
	print "Author  :\tLei Gao, UC,Riverside\n";
	print "Contact :\tLei Gao <highlei\@gmail.com>\n";
	print "\nUsage:	$0 [options]\n";
	print "\t-i	<str>	input the soap result file.";
	print " eg: blast_m8.blastn/soap.M0r2v0.soap\n";
	print "\t-p	<str>	input the position file.";
	print " eg: rice.genome.soap\n";
	print "\t-c	<str>	input the column.";
	print " [$col]\n";
	print "\t-f	<int>	input fore length.";
	print " [$fore]\n";
	print "\t-b	<int>	input back length.";
	print " [$back]\n";

	print "\n\t-h	display this help\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
	print "$0 -i blast_m8.blastn/soap.M0r2v0.soap -p genome.pos -c 1,2 -f $fore -b $back\n";
	print ("==========================| $0  end  |==================================\n\n");

    exit(0);
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
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}

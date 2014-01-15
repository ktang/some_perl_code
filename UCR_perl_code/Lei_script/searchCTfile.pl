#!/usr/bin/perl
# Copyright (c)  2009--
# Program:			searchCTfile
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.11.18
# Description:		extract seqs from ct file according to filelist or name
# version: 1.0s
# **************************
# Version: 1.0s
# **************************


# e-mail:highlei@itp.ac.cn
# ^-^

my $version="1.0s";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:s:l:r:t:c:");
my $infile		= $opt_i;
my $string		= $opt_s;
my $filelist	= (defined $opt_l) ? $opt_l : "";
my $reversive	= (defined $opt_r) ? $opt_r : "";
my $nameLen		= (defined $opt_t) ? $opt_t : 99999;
my $col		= (defined $opt_c) ? $opt_c-1 : 0;	## the column

if ($opt_h || $infile eq "" || ($string eq "" && $filelist eq "")) {
	usage();
}
use FileHandle;
use strict;

my $head = 100;
my $tail = 100;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$flag,$line,$key,$end);
my (%gnm);
my $kk	= 0;
my $skip	= "";


#===========================================================================================================
#====================                  main
#===========================================================================================================

if ($string ne "") {
	if (exists($gnm{$string})) {
		print STDERR "the string occupy many time!!\n";
	} else {
#		$string =~ s/(\W)/\\$1/g;
		$gnm{$string} = 1;	$kk++;
	}
}
for ($k1 = 0;$k1 < $col ;$k1++) {
	$skip	.= "\\S+\\s+";
}
print STDERR "skip=$skip\n";
if ($filelist ne "") {
	$file = new FileHandle("$filelist") || die("Cannot open listfile $filelist!\n");
	while (<$file>) {
		$_ =~s/^[\s|\t]+//g;	#print STDERR $_,",";
#		$_ =~s/[\s|\t|\r|\n]+$//g;
		$_=~/^$skip(\S+)/;		#print STDERR $1,"\n";
		if (exists($gnm{$1})) {
			print STDERR "the name:$_ occur many times!!\n";
		} else {
#			$_ =~ s/(\W)/\\$1/g;
			$gnm{$1} = 1;	$kk++;
		}
	}
	close $file || die;
}

print STDERR "\nread list is completed:$kk in file $filelist!\n";
#exit(0);
$m = 0;	$flag	= 0;
$file = new FileHandle ("$infile") || die("Cannot open file $infile!");
$i	= <$file>;
while (!eof) {
#for (; ;) {
	$i	=~/\=\s+\S+\s+(\S+)\s/;
#	foreach $key (keys(%gnm)) {
#		$key	= quotemeta($key);
		if (exists($gnm{$1}) && $reversive eq "") {
			if (length($i) > $nameLen) {
				print substr($i,0,$nameLen);
			} else {
				print $i;
			}
			$m++;	$flag	= 2;
			while (<$file>) {
				if ($_=~/\=/) {
					$flag	= 1;
					$i	= $_;
					last;
				}
				print $_;
			}
		} elsif (exists($gnm{$1})) {
			$flag	= 2;
			while (<$file>) {
				if ($_=~/\=/) {
					$flag	= 1;
					$i	= $_;
					last;
				}
			}
		}
#	}
	if ($flag	== 0 && $reversive ne "") {
		if (length($i) > $nameLen) {
				print substr($i,0,$nameLen);
			} else {
				print $i;
			}
		$m++;	$flag	= 2;
		while (<$file>) {
			if ($_=~/\=/) {
				$flag	= 1;
				$i	= $_;
				last;
			}
			print $_;
		}
		$m%100000!=0 || print STDERR $m," ";
	} elsif ($flag	== 0) {
		while (<$file>) {
			if ($_=~/\=/) {
				$i	= $_;
				last;
			}
		}
	} 
	if ($flag	== 2) {
		last;
	} else {
		$flag	= 0;
	}
}

close $file || die("Wrong!");

print STDERR "\nthere is $m sequences in this file\n";

print STDERR ("============================================================\n");
my $Time_End = sub_format_datetime(localtime(time()));
print STDERR "Running from [$Time_Start] to [$Time_End]\n";
$end = time();
printf STDERR ("Total execute time : %.2f s\n",$end-$start);
print STDERR ("............................................................\n");

#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################


sub usage
{
	print "Contact : Gaolei <highlei\@itp.ac.cn>";
	print "\nProgram : $0\nVersion: $version\n";
	print "Usage:\n	$0 \n";
	print "-i	input database file\n";
	print "			eg: syn.faa\n";
	print "-s	input querry sting\n";
	print "			eg: ARGTK\n";
	print "-l	input sequence name list file\n";
	print "			eg: seqname.lis\n";
	print "-r	reverse the result!\n";
	print "			eg: -r n; default: do not reverse\n";
	print "-t	the maxmum length of annotation line (name line)!\n";
	print "			eg: -t 1200; default: 999\n";
	print "-c	use the nth column in -l filelist\n";
	print "			eg: -c 2; default: 1, ie the 1st column\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i syn.faa -s AHKJL -l seqname.lis -r n -t 999\n";
	print "$0 -i syn.faa -l seqname.lis -c 2 -r n -t 999\n";
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



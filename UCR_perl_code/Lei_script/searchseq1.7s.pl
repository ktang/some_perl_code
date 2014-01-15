#!/usr/bin/perl
# Copyright (c)  2003--
# Program:			searchseq
# Author:			Gaolei <highlei@itp.ac.cn>
# Program Date:		2003.12.12
# Modifier:			Gaolei <highlei@itp.ac.cn>
# Last Modified:	2010.08.25
# Description:		extract seqs from database file according to filelist or name
# **************************
# Version: 1.1	add parameter -r
# Version: 1.2	match_quotemeta(),deal with special letter.
# Version: 1.2s	special version for only string (list)=sequence name.
# Version: 1.3s new parameter $nameLen, limit the length of annotation line (name line).
# Version: 1.4s	add the new parameter $col
# Version: 1.5s	use the fasta file as the filelist
# Version: 1.6s	removing the same name.
# Version: 1.7s	use the list/fasta description
# **************************


# e-mail:highlei@itp.ac.cn
# ^-^

my $version="1.7s";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:s:l:f:r:t:c:n:d:");
my $infile		= $opt_i;
my $string		= $opt_s;
my $filelist	= (defined $opt_l) ? $opt_l : "";
my $fastalist	= (defined $opt_f) ? $opt_f : "";
my $reversive	= (defined $opt_r) ? $opt_r : "";
my $nameLen		= (defined $opt_t) ? $opt_t : 99999;
my $col			= (defined $opt_c) ? $opt_c-1 : 0;	## the column
my $sameName	= (defined $opt_n) ? $opt_n : 0;	# 0: no removing; 1: removing
my $descript	= (defined $opt_d) ? $opt_d : 0;	# 0: use the $infile's; 1 use the fasta/list's

if ($opt_h || $infile eq "" || ($string eq "" && $filelist eq "" && $fastalist eq "" )) {
	usage();
}
use FileHandle;
use strict;

my $head = 100;
my $tail = 100;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$flag,$line,$key,$end);
my (%gnm,%name);
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
		if ($_=~/^$skip(\S+)\s+(.+)$/) {
			$k1	= $1;	$k2	= $2;
		} elsif ($_=~/^$skip(\S+)/) {
			$k1	= $1;	$k2	= " ";
		}
		if (exists($gnm{$k1})) {
			print STDERR "the name:$_ occur many times!!\n";
		} else {
#			$_ =~ s/(\W)/\\$1/g;
			$gnm{$k1} = $k2;	$kk++;
		}
	}
	close $file || die;
	print STDERR "\nread list is completed:$kk in file $filelist!\n";
}
if ($fastalist ne "") {
	$file = new FileHandle("$fastalist") || die("Cannot open fastalist $fastalist!\n");
	while (<$file>) {
		$_ =~s/^[\s|\t]+//g;	#print STDERR $_,",";
		if ($_=~/^>/) {
					#print STDERR $1,"\n";
			if ($_=~/^>$skip(\S+)\s+(.+)$/) {
				$k1	= $1;	$k2	= $2;
			} elsif ($_=~/^>$skip(\S+)/) {
				$k1	= $1;	$k2	= " ";
			}
			if (exists($gnm{$k1})) {
				print STDERR "the name:$_ occur many times!!\n";
			} else {
	#			$_ =~ s/(\W)/\\$1/g;
				$gnm{$k1} = $k2;	$kk++;
			}
		}
	}
	close $file || die;
	print STDERR "\nread fastalist is completed:$kk in file $fastalist!\n";
}
#exit(0);
$m = 0;	$flag	= 0;
$file = new FileHandle ("$infile") || die("Cannot open file $infile!");
$i	= <$file>;
while (!eof) {
#for (; ;) {
	$i	=~/>(\S+)\s/;
#	foreach $key (keys(%gnm)) {
#		$key	= quotemeta($key);
		if (exists($gnm{$1}) && $reversive eq "") {
			if ($sameName == 1 && exists($name{$1})) {
				while (<$file>) {
					if ($_=~/^\>/) {
						$flag	= 1;
						$i	= $_;
						last;
					}
				}
				next;
			} else {
				$name{$1}	= 1;
			}
			if (length($i) > $nameLen) {
				print substr($i,0,$nameLen);
			} else {
				if ($descript != 0) {
					$i	=~/>(\S+)\s/;
					print ">$1 ",$gnm{$1},"\n";
				} else {
					print $i;
				}
			}
			$m++;	$flag	= 2;
			while (<$file>) {
				if ($_=~/^\>/) {
					$flag	= 1;
					$i	= $_;
					last;
				}
				print $_;
			}
		} elsif (exists($gnm{$1})) {
			$flag	= 2;
			while (<$file>) {
				if ($_=~/^\>/) {
					$flag	= 1;
					$i	= $_;
					last;
				}
			}
		}
#	}
	if ($flag	== 0 && $reversive ne "") {
		if ($sameName == 1 && exists($name{$1})) {
				while (<$file>) {
					if ($_=~/^\>/) {
						$flag	= 1;
						$i	= $_;
						last;
					}
				}
				next;
		} else {
				$name{$1}	= 1;
		}
		if (length($i) > $nameLen) {
				print substr($i,0,$nameLen);
			} else {
				if ($descript != 0) {
					$i	=~/>(\S+)\s/;
					print ">$1 ",$gnm{$1},"\n";
				} else {
					print $i;
				}
			}
		$m++;	$flag	= 2;
		while (<$file>) {
			if ($_=~/^\>/) {
				$flag	= 1;
				$i	= $_;
				last;
			}
			print $_;
		}
		$m%100000!=0 || print STDERR $m," ";
	} elsif ($flag	== 0) {
		while (<$file>) {
			if ($_=~/^\>/) {
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
	print "-f	input sequence name list file fasta format\n";
	print "			eg: seqname.fasta\n";
	print "-r	reverse the result!\n";
	print "			eg: -r n; default: do not reverse\n";
	print "-t	the maxmum length of annotation line (name line)!\n";
	print "			eg: -t 1200; default: 999\n";
	print "-c	use the nth column in -l filelist\n";
	print "			eg: -c 2; default: 1, ie the 1st column\n";
	print "-n	removing the same name, output unique sequences.\n";
	print "			eg: -n 1; default: 0; 0: no removing; 1: removing\n";
	print "-d	use the infile or fasta description.\n";
	print "			eg: -d 1; default: 0; 0: use $infile's; 1: use fasta's\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i syn.faa -s AHKJL -l seqname.lis -r n -t 999\n";
	print "$0 -i syn.faa -l seqname.lis -c 2 -r n -t 999\n";
	print"\n============================================================\n\n";

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



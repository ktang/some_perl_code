#!/usr/bin/perl
# Copyright (c) Lei Gao 2009
# Program:			blast2histogram
# Author:			Gaolei <highlei@gmail.com>
# Program Date:		2009.08.24
# Modifier:			Gaolei <highlei@gmail.com>
# Last Modified:	2009.08.24
# Description:	read blast result file(m=0) to get the data for histogram. $len = 0:whole length match; -1: 1 base shorter than whole length
#				(hit length or $overlap_len);
#**************************
# Version: 1.1	output the sequence of hit
#**************************
#refer to blastTable1.1.pl

# e-mail:highlei@gmail.com

my $version="1.1";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:o:l:d:q:e:s:t:p:c:y:");
my $infile	= (defined $opt_i) ? $opt_i : "";	#"pn.BLASTP";  #结果文件
my $outfile	= (defined $opt_o) ? $opt_o : "";
my $database	= (defined $opt_d) ? $opt_d : "";
my $query	= (defined $opt_q) ? $opt_q : "";
my $evalue	= (defined $opt_e) ? $opt_e : 10;          #BlAST 搜索的域值
#my $lisfile	= $opt_l;
my $len		= (defined $opt_l) ? $opt_l : 0;	# 0:whole length match; -1: 1 base shorter than whole length
my $hsp_1st	= (defined $opt_s) ? $opt_s : 1;	# default: output 1st hsp. 0: output all hsp
my $hit_1st	= (defined $opt_t) ? $opt_t : -1;	# default: output all hits.
my $overlap_len	= (defined $opt_p) ? $opt_p : 100;
my $score	= (defined $opt_c) ? $opt_c : 0;
my $identity	= (defined $opt_y) ? $opt_y : 100;

if ($opt_h || $infile eq "") {# || ($from eq "" xor $to eq "")){
	usage();
}
use Bio::SearchIO;
use FileHandle;
use strict;

my $head = 100;
my $tail = 100;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$line,$count,$flag,$block,$a,$b,$end);
my (@buf,@tmp,@buf1,@buf2,@op);
my $ii=0;
my $jj=0;
my $kk=0;
my %entry=();
my %query=();
my %gnm=();
my $key="";
my ($r,$h,$hsp,$qlen,$dblen);
my ($qname, $dbname, $qidentity, $alnLen, $qstart, $qend, $dbstart, $dbend, $qevalue, $qbits);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#print STDERR "\-l=$len\n";
$m	= 0;
############################################ read genome files ######################################################
if ($query ne "") {
	$i	= -1;
	$file = new FileHandle ("$query") || die("Cannot open the file $query!\n");
	while (<$file>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		if ($_ =~/^>(\S+)/) {
		#	if (exists($genome{$1})) {
		#		die("There are something wrong in $gnmFiles:$1\n");
		#	} else {
				$i	= $1;
				$query{$i}	= "";
				$m++;
		#	}
		} else {
			$query{$i} .= $_;
		}
	}
	close $file || die("Wrong!");
	print STDERR "There are $m sequence in file $query\n";
}
$m	= 0;
############################################ read genome files ######################################################
if ($database ne "") {
	$i	= -1;
	$file = new FileHandle ("$database") || die("Cannot open the DB $database!\n");
	while (<$file>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		if ($_ =~/^>(\S+)/) {
		#	if (exists($genome{$1})) {
		#		die("There are something wrong in $gnmFiles:$1\n");
		#	} else {
				$i	= $1;
				$gnm{$i}	= "";
				$m++;
		#	}
		} else {
			$gnm{$i} .= $_;
		}
	}
	close $file || die("Wrong!");
	print STDERR "There are $m sequence in file $database\n";
}

$m = 0;	$n	= 0;
if ($outfile ne "") {
	open(OUT, ">$outfile") || die("Can not open file: $outfile\n");
}
$k1	= 0;	$flag	= 0;	$a = 0;
$file = new FileHandle ("$infile") || die("Cannot open file: $infile");

while(<$file>)
{
	if ($_ =~/^(\S+)\s+(\S+)\s+([\d|\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([\d|\.]+)/) {
		last if $9 >= $evalue;	$b	= $_;
		if (!exists($entry{$1})) {
			if ($a != 0) {
				print ">",$qname,"\n";
				for ($k2 = 0;$k2 < $qlen ;$k2++) {
					print $buf[$k2]," ";
				}
				print "\n";
			} else {
#				print STDERR "a=$a---\n";
			}

			$entry{$1}{"num"}	= 1;#	$flag	= 0;
			$qlen	= length($query{$1});
			$dblen	= length($gnm{$2});	splice(@buf,0);
			for ($k2 = 0; $k2 < $qlen ;$k2++) {
				$buf[$k2]	= 0;
			}
			$kk	= 0;	$k1	= 0;	$a	= 0;	$ii++;
		} else {
			$dblen	= length($gnm{$2});
			if ($hit_1st != -1 && $kk > $hit_1st) { 
#				print STDERR "kk=$kk---\n";
				next;
			}
		}
		$qname=$1; $dbname=$2; $qidentity=$3; $alnLen=$4; $qevalue=$9; $qbits=$10;
		$dbstart	= $7 > $8 ? $8 : $7;	$dbend		= $7 > $8 ? $7 : $8;
		$qstart		= $5 > $6 ? $6 : $5;	$qend		= $5 > $6 ? $5 : $6;
#		$qstart=$5; $qend=$6; $dbstart=$7; $dbend=$8; 
		
		$dbname	=~ /\_(\d+)x?$/;
		$k3	= $1;
		if ($qevalue	== 0 && $qbits < 50) {
			$kk--;	#print STDERR "qevalue=$qevalue---\n";
			next;
		}
		if (exists($entry{$qname}{$dbname})) {
			$jj++;
			$entry{$qname}{$dbname}++;
			if ($hsp_1st == 0) {
			} elsif ($entry{$qname}{$dbname} > $hsp_1st){ #$hsp_1st == 1 && $jj > 1) {
#				print STDERR "jj=$jj---\n";
				next;
			}
		} else {
			$jj	= 0;
			$entry{$qname}{$dbname}	= 0;
		}
		if ($qidentity < $identity) {
#			print STDERR "qidentity=$qidentity---\n";
			next;
		}
		if ($len == 0) {
			if ($alnLen >= $overlap_len || $alnLen == $dblen) { # hsp length == hit length (database)
			} else {
#				print STDERR "alnLen=$alnLen,dblen=$dblen---\n";
				next;
			}
		} elsif ($len < 0) {
			if ($alnLen >= $dblen+$len) {
			} else {
#				print STDERR "alnLen=$alnLen-<0\n";
				next;
			}
		} elsif ($alnLen < $len) {
#			print STDERR "alnLen=$alnLen-<\n";
			next;
		}
		if ($kk == 0) {
			$n++;	$kk++;
		} else {
			$kk++;
		}
		for ($k2 = $qstart -1; $k2 < $qend ;$k2++) {
			$buf[$k2]	+= $k3;	$a	+= $k3;
		}
		if ($jj == 1 && $outfile ne "") {
			print OUT $qname,"\t",$dbname,"\n";
#			if ($database ne "") {
#				if (exists($gnm{$dbname})) {
#					print OUT $gnm{$dbname};
#				} else {
#					print STDERR "no seq:",$dbname,"\t",$qname,"\n";
#				}
#			}
		}
	$m++;	#	print STDERR $b;
	} else {
		print STDERR "the format is wrong:",$_;
	}
}
if ($a != 0) {
	print ">",$qname,"\n";
	for ($k2 = 0;$k2 < $qlen ;$k2++) {
		print $buf[$k2]," ";
	}
	print "\n";
}

close $file || die;

#print	STDERR "\nread file $infile completely,element num=$k1!\n";


if ($outfile ne "") {
	close(OUT) || die("Can not close file: $outfile\n");
}
#--------------------------------------------------------------------------------------------------------
#print STDERR "\n\nBio::SearchIO end!\n";
print STDERR "\n\nthere are $ii,$n,$m allquery,hitquery,pairs in file: $infile\n";
#---------------------------------------output-----------------------------------------------




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
	print "-i	input file\n";
	print "			eg: amo.blast\n";
	print "-o	output file including sequences found in genomes\n";
	print "			eg: -o align_seq.fna\n";
	print "-e	evalue for blast significance (or evalue)\n";
	print "			eg: 0.01;default: 10 \n";  
	print "-l	cutoff of length of alignment lengh. 0:whole length match; minus: i.e. -1: 1 base shorter than whole length\n";
	print "			eg: -l 30; default: 0 \n"; 
	print "-p	the overlap length of neighbor segments\n";
	print "			eg: -p 100;default: 100\n";
	print "-s	control about whether only output the 1st hsp or multi\n";
	print "			eg: -s 2; default: -s 1,i.e. only output the 1st hsp \n";
	print "-t	control about whether only output the 1st hit or multi\n";
	print "			eg: -t 2; default: -t -1,i.e. output all hit \n";
	print "-y	cutoff of identity\n";
	print "			eg: -y 99; default: 0\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i amo.tblastn -t 1000 > blastTable.txt\n";
	print "$0 -i amo.tblastn -e 10 -t -1 -l -1 -y 100 > blastTable.t1y00l-1.txt\n";
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
	exit(0);

}


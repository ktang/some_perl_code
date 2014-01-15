#!/usr/bin/perl
# Copyright (c)  2008-
# Program:			qsub4miRNAInGnm
# Author:			Gaolei <highlei@gmail.com>
# Program Date:		2008.11.13
# Modifier:			Gaolei <highlei@gmail.com>
# Last Modified:	2010.03.03
# Description:  qsub jobs.
#**************************
# Version: 1.1	add parameter $op_P_F; not read all query into memory, yield one sh file when read part of the sequences.
#**************************


# e-mail:highlei@gmail.com

my $version="1.1";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hp:s:b:q:o:m:d:i:r:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $qsub		= (defined $opt_b) ? $opt_b : "";#qsub";
my $pro			= (defined $opt_p) ? $opt_p : "blastall -p blastx";
my $database	= $opt_d;
my $query		= $opt_q;
my $param		= $opt_m;
my $splitQuery	= (defined $opt_s) ? $opt_s : 0;
my $queryparam	= (defined $opt_i) ? $opt_i : "-i";
my $op_P_F		= (defined $opt_o) ? $opt_o : "-o " . "__qsub_$start";
my $separator	= (defined $opt_r) ? $opt_r : ">";	#"BLASTN";

if ($opt_h || $pro eq ""){# || $database eq "" || $query eq "") {
	usage();
}
use FileHandle;
use strict;



my ($i,$j,$k,$m,$n,$sh,$k1,$k2,$k3,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@bufsh,@tmpop,@output);
my $key="";
my $flag=0;
my ($bfile,$sfile);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\nbatch version 1.3\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n b  %35s : %-25s","qsub or not",$qsub);#%35s
	printf STDERR ("\n p  %35s : %-25s","the program name",$pro);#%35s
	printf STDERR ("\n d  %35s : %-25s","database",$database);
	printf STDERR ("\n q  %35s : %-25s","query",$query);
	printf STDERR ("\n o  %35s : %-25s","output parameter+output files",$op_P_F);
	printf STDERR ("\n i  %35s : %-25s","query parameters",$queryparam);
	printf STDERR ("\n m  %35s : %-25s","other parameters",$param);
	printf STDERR ("\n s  %35s : %-25s","split query",$splitQuery);
	printf STDERR ("\n r  %35s : %-25s","separator to split the query",$separator);
	printf STDERR ("\n x  %35s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "b") {print STDERR "please input qsub or not\n"; $qsub	= <STDIN>;	$qsub	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "p") {print STDERR "please input the program name:\n"; $pro	= <STDIN>;	$pro	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "d") {print STDERR "please input the database:\n";$database	= <STDIN>;$database	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "q") {print STDERR "please input the queryeters:\n";$query	= <STDIN>;$query	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "o") {print STDERR "please input the output parameter+output files:\n";$op_P_F	= <STDIN>;$op_P_F	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "i") {print STDERR "please input the query parameters:\n";$queryparam	= <STDIN>;$queryparam	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "m") {print STDERR "please input other parameters:\n";$param	= <STDIN>;$param	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input how many files split query to:\n";$splitQuery	= <STDIN>;$splitQuery	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "r") {print STDERR "please input the separator to split the query: e.g.:fasta:\">\",blast:\"BLAST\",each line:\"\"\n";$separator	= <STDIN>;$separator	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(1);}
}

#$param	=~s/&+/ /g;
$database	=~s/\"$//g;	$database	=~s/^\"//g;
$query	=~s/\"$//g;	$query	=~s/^\"//g;
$op_P_F	=~s/\"$//g;		$op_P_F	=~s/^\"//g;
#$queryparam	=~s/\"$//g;	$queryparam	=~s/^\"//g;
$param	=~s/\"$//g;		$param	=~s/^\"//g;

#print STDERR "\npro=$pro\tparam=$param\tbachfile=$database\top=$output\tdir=$dir\n";

############################################ read file ######################################################
$i	= 0;	$k3	= 0;	$m	= 0;	$sh	= 0;
if ($query ne "") {
	if ($query =~/^(.+\/)[^\/]*\*([^\*]+)$/) {
		$k	= $2;	$j	= $1;
	} elsif ($query =~/^[^\/]*.*\*([^\*]+)$/) {
		$k	= $2;	$j	= "./";
	} elsif ($query =~/^(.+\/)([^\/]+)$/) {
		$k	= $2;	$j	= $1;
	} elsif ($query =~/^([^\*]+)$/) {
		$k	= $1;	$j	= "./";
	} else {
		print STDERR	$query;
		die("reinput $query!\n");
	}
	$k	= quotemeta($k);
	opendir(FDIR, $j) || die("Can not open dir: $j\n");
	while ($file=readdir(FDIR)) {
		if ($file=~/$k$/) {
			$bfile	= $j ."/$file";
			$i++;
			print STDERR "\n----------------------------------------------------------------------\nthe $i file:\t$bfile\n";
#--------------------------------------------------------------------------------------------------
if ($splitQuery > 1) {
	$sfile = new FileHandle ("$bfile") || die;
	$k1	= -1;	splice(@bufi,0);
	while(<$sfile>)
	{
		if ($separator eq "") {
			$k1++;
#			$bufi[$k1]	= $_;
		} elsif ($_ =~/^$separator/) {
			$k1++;
#			$bufi[$k1]	= $_;
		} else {
#			$bufi[$k1]	.= $_;
		}
	}
	close $sfile || die;
	print	STDERR	"\tsequence number = ",++$k1,"\n";

	$k2 = 1;
	$sfile = new FileHandle ("$bfile") || die;
	$bufo[$m]	= "__$file.$m.fa";
	open(SFILE,">$bufo[$m]") || die("Can not open file: $bufo[$m]\n");
	while (<$sfile>) {
		if ($separator eq "" || $_ =~/^$separator/) {
			if ($k2%($k1/$splitQuery) == 0) {
				my	$k10	= $_;
				if ($k1-$k2 < $splitQuery && $m > $splitQuery -2) {
					print SFILE $k10;
					next;
				}
				close(SFILE);

				$output[$k3]	= $op_P_F . ".$k3.output";
				print	STDERR "$k2\t$pro $param $database $queryparam $bufo[$m] $output[$k3]\n";
				$bufsh[$sh]	= "__$file.$sh.sh";
				open(SHF,">$bufsh[$sh]") || die("Can not open file: $bufsh[$sh]\n");
				print	SHF "$pro $param $database $queryparam $bufo[$m] $output[$k3]\n";
				close(SHF);
				if ($qsub eq "") {
					system("chmod 755 $bufsh[$sh]");
					system("./$bufsh[$sh] &");
				} elsif ($qsub eq "qsub") {
					system("$qsub $bufsh[$sh]");
				} else {
					system("$qsub $bufsh[$sh] &");
				}
				$n	= "";	$m++;	$sh++;	$k3++;	if($sh>101){die("$sh");}

				$bufo[$m]	= "__$file.$m.fa";
				open(SFILE,">$bufo[$m]") || die("Can not open file: $bufo[$m]\n");
				print SFILE $k10;
			} else {
				print SFILE $_;
			}
			$k2++;
#			$bufi[$k1]	= $_;
		} else {
			print SFILE $_;
		}
	}
	close(SFILE);
	$output[$k3]	= $op_P_F . ".$k3.output";
	print	STDERR "$k2\t$pro $param $database $queryparam $bufo[$m] $output[$k3]\n";
	$bufsh[$sh]	= "__$file.$sh.sh";
	open(SHF,">$bufsh[$sh]") || die("Can not open file: $bufsh[$sh]\n");
	print	SHF "$pro $param $database $queryparam $bufo[$m] $output[$k3]\n";
	close(SHF);
	if ($qsub eq "") {
		system("chmod 755 $bufsh[$sh]");
		system("./$bufsh[$sh] &");
	} elsif ($qsub eq "qsub") {
		system("$qsub $bufsh[$sh]");
	} else {
		system("$qsub $bufsh[$sh] &");
	}
	$n	= "";	$m++;	$sh++;	$k3++;	if($sh>101){die("$sh");}

	close $sfile || die;
} else {
	$output[$k3]	= $op_P_F . ".$k3.output";
	print	STDERR "\t$pro $param $database $queryparam $query $output[$k3]\n";
	$bufsh[$sh]	= "__$file.$sh.sh";
#	print	STDERR "__$file.$sh.sh\n";
	open(SHF,">$bufsh[$sh]") || die("Can not open file: $bufsh[$sh]\n");
	print	SHF "$pro $param $database $queryparam $query $output[$k3]\n";
	close(SHF);
	if ($qsub eq "") {
		system("chmod 755 $bufsh[$sh]");
		system("./$bufsh[$sh] &");
	} elsif ($qsub eq "qsub") {
		system("$qsub $bufsh[$sh]");
	} else {
		system("$qsub $bufsh[$sh] &");
	}
#	system("qsub $bufsh[$sh]");
	$sh++;	$k3++;	if($sh>101){die("$sh");}
}
#--------------------------------------------------------------------------------------------------
		}
	}
	closedir(FDIR);
}
print	STDERR "\n\nThere are $sh jobs which were submitted!\n\n";

for ($k1 = 0;$k1 < $m ;$k1++) {
#	system("rm $bufo[$k1]");
}
for ($k1 = 0;$k1 < $sh ;$k1++) {
#	system("rm $bufsh[$k1]");
}

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
	print "Contact : Gaolei <highlei\@gmail.com>";
	print "\nProgram : $0\nVersion: $version\n";
	print "Usage:\n	$0 \n";
	print "-b	input qsub or not\n";
	print "			eg: -b qsub; default: qsub\n";
	print "-p	input program name\n";
	print "			eg: \"blastall -p blastx\"\n";
	print "-d	input database\n";
	print "			eg: \"-d bacteria751+4plant.faa\"\n";
	print "-i	input query parameter\n";
	print "			eg: \"-i\"; default:\"-i\"\n";
	print "-q	input query files name\n";
	print "			eg: ../wgsLucyTrim_*.fna\n";
	print "-s	input how many files split the query into\n";
	print "			eg: -s 10; default:0,i.e. donnot split\n";
	print "-r	 input the separator to split the query: e.g.:fasta:\">\",blast:\"BLAST\"\n";
	print "			eg: -r \"\"; default:>\n";
    print "-m	input all other parameters.\n";
    print "			eg: \"-e 1e-5 -o wgs_15.blastx\" \n"; 
	print "-o	input the output parameter+output files:\n\n";
	print "			eg: -o  \"-o __qsub_rice_\"\n";
    print "-h	display this lines\n";
	print "		Note: please add quotation mark, if you input parameter in command line!\n";
    print "\nExample:\n";
    print "$0 -p \"blastall -p blastx\" -d bacteria+p.faa -i \"-i\" -q \"wgs_*.faa\" -m \"-e 1e-5 -o wgs_17.blastx\" -s 16\n";
    print "$0 -b \"\" -p \"blastall -p blastx\" -d bacteria+p.faa -i \"-i\" -q \"wgs_*.faa\" -m \"-e 1e-5 -o wgs_17.blastx\" -s 16\n";
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



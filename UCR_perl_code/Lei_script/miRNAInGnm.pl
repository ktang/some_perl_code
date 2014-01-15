#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.11.18
# Description:	read the blast m8 result and cut the segments from the genome files, then do unafold, then do ctStrucFilter
#**************************
# Version: 1.0
#**************************
# e-mail:highlei@gmail.com

my $version="1.0";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:g:l:p:o:q:c:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $gnmFiles	= (defined $opt_g) ? $opt_g : "";
my $len			= (defined $opt_l) ? $opt_l : "20,100,20,300";
my $path		= (defined $opt_p) ? $opt_p : "./";
my $output		= (defined $opt_o) ? $opt_o : 0;		# 0: output the .fa and .ct files which 
my $qsub_set	= (defined $opt_q) ? $opt_q : "-b \"\" -s 2" ; # set -s and -r
my $ctStrucFlt	= (defined $opt_c) ? $opt_c : "-m 4 -n 2 -v 0.5";

if ($opt_h){
	usage();
}

sub numerically{$a<=>$b};

use FileHandle;
use strict;



my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$k4,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@genome,@gnmName,@gnmLen);
my (%gnm,%seg);
my $key="";
my ($endLen,$minLen,$addLen,$maxLen,$foldFileNum);
my ($foldFile,$foldFileOut);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %45s : %-25s","input blast result (m=8) file",$infile);#%45s
	printf STDERR ("\n g  %45s : %-25s","input genome file",$gnmFiles);
	printf STDERR ("\n l  %45s : %-25s","input length",$len);
	printf STDERR ("\n p  %45s : %-25s","input the path of perl scripts",$path);#%45s
	printf STDERR ("\n o  %45s : %-25s","output format",$output);
	printf STDERR ("\n q  %45s : %-25s","the parameter set for qsub_noFasta.pl",$qsub_set);
	printf STDERR ("\n c  %45s : %-25s","the parameter set for ctStructureFilter.pl",$ctStrucFlt);
#	if($zero==1) {printf STDERR ("\n z  %45s : %-25s","output coverage region?","1");}
#	elsif ($zero == 0) {printf STDERR ("\n z  %45s : %-25s","concise output","0");}
#	else {printf STDERR ("\n z  %45s : %-25s","output all genome","2");}
	printf STDERR ("\n x  %45s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n------------------------------------------------------------\n\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input the blast (m=8) result file:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "g") {print STDERR "please input the genome file:\n"; $gnmFiles	= <STDIN>;$gnmFiles	=~s/[\s|\t|\r|\n]+$//g;}
#	elsif($yesorno eq "z") {$zero		= ($zero+1)%3;}
	elsif($yesorno eq "l") {print STDERR "please input len(endLen,minLen,addLen,maxLen):\n";$len	= <STDIN>;$len	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "p") {print STDERR "please input the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl):\n";$path	= <STDIN>;$path	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "o") {print STDERR "please input output format:\n";$output	= <STDIN>;$output	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "q") {print STDERR "please input the parameter set for qsub_noFasta.pl:\n";$qsub_set	= <STDIN>;$qsub_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "c") {print STDERR "please input the parameter set for ctStructureFilter.pl:\n";$ctStrucFlt	= <STDIN>;$ctStrucFlt	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(0);}
}

if ($len =~ /(\d+)\,(\d+)\,(\d+),(\d+)/) {
	$endLen	= $1;	$minLen	= $2;
	$addLen	= $3;	$maxLen	= $4;
	print STDERR "endLen=$endLen\tminLen=$minLen\taddLen=$addLen\tmaxLen=$maxLen\n";
} else {
	print STDERR "please input the correct len!\n";
	usage();
}
#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";

############################################ read file ######################################################


$k1	= 0;
$file = new FileHandle ("$infile") || die("Cannot open file: $infile");

while(<$file>)
{
	if ($_ =~/^^(\S+)\s+(\S+)\s+([\d|\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+/) {
		$j	= $2;
		$m	= $7 > $8 ? $8 : $7;
		$n	= $7 > $8 ? $7 : $8;
		if (exists($gnm{$j}{$m}{$n})) {
			$gnm{$j}{$m}{$n}	.= ",".$1;
#			print STDERR "This occurs more than one time:",$_;
		} else {
			$gnm{$j}{$m}{$n}	= $1;
		}
		$k1++;
	} else {
		print STDERR "the format is wrong:",$_;
	}
}

close $file || die;

print	STDERR "\nread file $infile completely,element num=$k1!\n";
############################################ read genome files ######################################################
if ($gnmFiles ne "") {
	$i	= -1;
	$file = new FileHandle ("$gnmFiles") || die;
	while (<$file>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		if ($_ =~/^>(\S+)/) {
		#	if (exists($genome{$1})) {
		#		die("There are something wrong in $gnmFiles:$1\n");
		#	} else {
				$i++;
				$genome[$i]	= "";
				$gnmName[$i]	= $1;	#$k1++;
		#	}
		} else {
			$genome[$i] .= $_;
		}
	}
	close $file || die("Wrong!");
	print STDERR "\nThere are ",$i+1," sequence in file $gnmFiles\n";
	############################################ chromosome legth ######################################################
	$m	= 1;
	for ($k2=0; $k2 <= $i;$k2++) {
		$gnmLen[$k2]	= length($genome[$k2]);
		print STDERR "the ",$m,"th sequence: length = $gnmLen[$k2] bp\n";
		$m++;
	}
	print STDERR "\n";
}

############################################ creat %seg and output ######################################################

$foldFileNum	= 0;
system("rm $infile.miRNAInGnm*");

$k3	= 0;
for ($k1 = 0;$k1 <= $i ;$k1++) {
	$j	= $gnmName[$k1];
	foreach $m (sort numerically keys %{$gnm{$j}} ) {				# m	= start
		$k4 = 0;	$line	= "";
		foreach $n (sort numerically keys %{$gnm{$j}{$m}} ) {
			$k4 = $k4 > $n ? $k4 : $n;
			$line	.= $gnm{$j}{$m}{$n}.",";
		}
		$n	= $k4;
#		foreach $n (sort numerically keys %{$gnm{$j}{$m}} ) {		# n	= end
			$a = $m > $endLen ? $m-$endLen : 0;
			$b	= $a + $maxLen > $gnmLen[$k1]? $gnmLen[$k1]-$a : $maxLen;
			for ($k2 = $minLen; $k2 <= $b ;$k2+=$addLen) {
				$key	= "$j\_$a\_$k2";
				if (exists($seg{$key})) {
					$seg{$key}->[0]	.= "\t".$line.";$m,$n;".($m-$a).",".($n-$a);
				} else {
					$seg{$key}->[0]	= $line.";$m,$n;".($m-$a).",".($n-$a);
					$seg{$key}->[1]	= substr($genome[$k1],$a,$k2);	$k3++;
				}
	#			print ">$j\_$a\_$k2\t",$gnm{$j}{$m}{$n},"\n",substr($genome[$k1],$a,$k2),"\n";
			}
			$a	= $n + $endLen > $gnmLen[$k1]? $gnmLen[$k1] : $n+$endLen;
			$b	= $a - $maxLen > 0 ? $maxLen : $a;
			for ($k2 = $minLen;$k2 <= $b ;$k2+=$addLen) {
				$key	= "$j\_" . ($a-$k2) . "\_$k2";
				if (exists($seg{$key})) {
					$seg{$key}->[0]	.= "\t".$line.";$m,$n;".($m-$a+$k2).",".($n-$a+$k2);
				} else {
					$seg{$key}->[0]	= $line.";$m,$n;".($m-$a+$k2).",".($n-$a+$k2);
					$seg{$key}->[1]	= substr($genome[$k1],$a-$k2,$k2);	$k3++;
				}
	#			print ">$j\_",$a-$k2,"\_$k2\t",$gnm{$j}{$m}{$n},"\n",substr($genome[$k1],$a-$k2,$k2),"\n";
			}
#		}
	}
	$k4	= 0;
	$foldFile		= "__Tunafold_hybrid-ss-min_sDAT_$start.fa";		# $foldFileNum
	$foldFileOut	= "__Tunafold_hybrid-ss-min_sDAT_$start.fa.out";
	open(FOLDTMPF, ">$foldFile") || die("Can not open file: $foldFile\n");
	foreach $key (sort keys %seg) {
		print FOLDTMPF ">$key\t",$seg{$key}->[0],"\n",$seg{$key}->[1],"\n";	$k4++;
	}
	close(FOLDTMPF);
	print STDERR "\nthere are $k4 segments in $j -------------------------------------------------------\n\n";
#	if ($j =~/[5|0]/) {
	print STDERR "\tperl $path/qsub_nonFasta1.1.pl"."\n";
	system("perl $path/qsub_nonFasta1.1.pl -0 0 $qsub_set -p hybrid-ss-min -i \"\" -m \"-s DAT\" -q $foldFile -o \"> $foldFileOut\" 2> $foldFile.err");	# --mfold=*,*,5 
	for ($k2 = 0;$k2 < $k4 ;) {
		sleep(30);
		system("less __$foldFile*ct | gawk '/=/' | wc > $foldFile.tmp");
		open(NUMTMPF, "$foldFile.tmp") || die("Can not open file: $foldFile\n");
		$a	= <NUMTMPF>;
		close(NUMTMPF);
#		print STDERR "\twc=$a";
		$a	=~ /(\d+)\s+(\d+)\s+(\d+)/;
		$k2	= $1;
#		print STDERR "\tunafold: $k4\t$k2\n";
	}
	print STDERR "\tcat __$foldFile*ct > $foldFile.ct\n";
	system("cat __$foldFile*ct > $foldFile.ct");
	print STDERR "\tcat $foldFileOut* > $foldFileOut\n";
	system("cat $foldFileOut*put > $foldFileOut");
	
	print STDERR "\tperl $path/ctStructureFilter4miRNAInGnm.pl \n";
	system("perl $path/ctStructureFilter4miRNAInGnm.pl -0 0 $ctStrucFlt -i $foldFile.ct -o $foldFile.stem > $foldFile.ctFlt2.2.txt 2> $foldFile.err");
	print STDERR "\tperl $path/searchseq1.4s.pl\n";
	system("perl $path/searchseq1.4s.pl -i $foldFile -l $foldFile.ctFlt2.2.txt > $foldFile.ctFlt2.2.fa 2> $foldFile.err");
	print STDERR "\tperl $path/searchCTfile.pl\n";
	system("perl $path/searchCTfile.pl -i $foldFile.ct -l $foldFile.ctFlt2.2.txt > $foldFile.ctFlt2.2.fa.ct 2> $foldFile.err");
	system("cat $foldFile.ctFlt2.2.fa >> $infile.miRNAInGnm.fa");
	system("cat $foldFile.ctFlt2.2.fa.ct >> $infile.miRNAInGnm.fa.ct");
	system("cat $foldFile.ctFlt2.2.txt >> $infile.miRNAInGnm.txt");
	system("rm __$foldFile* $foldFile*"); #$foldFileOut*put $foldFile.tmp
#		system("perl ");
#	} else {
#		system("hybrid-ss-min -s DAT $foldFile > $foldFileOut");
#	}
#	$foldFileNum++;
	%seg=();
#	sub_end_program();
}

print STDERR "there are $k3 segments\n";
############################################ output ######################################################

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
	print "-g	input the genome file\n";
	print "			eg: genome.fa\n";
	print "-l	please input len(endLen,minLen,addLen,maxLen)\n";
	print "			eg: -l 10,75,50,500; default: $len\n";
	print "-p	the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl)\n";
	print "			eg: -p ../soft; default: $path, i.e. current directory\n";
	print "-o	output format\n";
	print "			eg: -o 1; default: $output\n";
	print "-q	the parameter set for qsub_noFasta.pl\n";
	print "			eg: -q \"\"; default: $qsub_set\n";
	print "-c	the parameter set for ctStructureFilter.pl\n";
	print "			eg: -c \"\"; default: $ctStrucFlt\n";
	print "-h	display this lines\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
#	print "$0 -i blast_m8.blastn -t 90\n";
	print "$0 -i blast_m8.blastn -g genome.fa -l $len -p $path -o $output -q $qsub_set -c $ctStrucFlt\n";
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
	return;
}

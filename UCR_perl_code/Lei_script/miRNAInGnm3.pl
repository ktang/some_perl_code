#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.11.18
# Description:	read the blast m8 result and cut the segments from the genome files, then do unafold, then do ctStrucFilter
#**************************
# Version: 2.0	do ctFilter for each qsub results. use do_fold_ctFlt4miRNAInGnm.pl
#**************************
# e-mail:highlei@gmail.com

my $version="3.0";
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
my $fold_ctFlt	= (defined $opt_c) ? $opt_c : "-m 4 -n 2";

if ($opt_h){
	usage();
}

my $qsub		= "qsub4miRNAInGnm.pl -0 0";
my $ctFlt		= "do_fold_ctFlt4miRNAInGnm.pl -0 0";#4miRNAInGnm
my $searchSeq	= "searchseq1.4s.pl";
my $searchCT	= "searchCTfile.pl";

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
	printf STDERR ("\n c  %45s : %-25s","the parameter set for do_fold_ctFlt4miRNAInGnm.pl",$fold_ctFlt);
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
	elsif($yesorno eq "c") {print STDERR "please input the parameter set for do_fold_ctFlt4miRNAInGnm.pl:\n";$fold_ctFlt	= <STDIN>;$fold_ctFlt	=~s/[\s|\t|\r|\n]+$//g;}
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
	if ($_ =~/^(\S+)\s+(\S+)\s+([\d|\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+/) {
		$j	= $2;
		$m	= $7 > $8 ? $8 : $7;
		$n	= $7 > $8 ? $7 : $8;
		if (exists($gnm{$j}{$m}{$n})) {
			$i	= $1;
			$i	=~ /\_(\d+)x*$/;	$k4	= $1;
			$gnm{$j}{$m}{$n}	=~ /\_(\d+)x*$/;	$k3	= $1;	#	print STDERR "old=$k3,new=$k4\n";
			if ($k4 < $k3) {
				$gnm{$j}{$m}{$n}	= $i.",".$gnm{$j}{$m}{$n};
			} else {
				$gnm{$j}{$m}{$n}	.= ",".$i;
			}													#	print STDERR $gnm{$j}{$m}{$n},"\n";
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
system("rm $infile.miRiG*");

$k3	= 0;	my $k6	= 0;
for ($k1 = 0;$k1 <= $i ;$k1++) {
	$j	= $gnmName[$k1];	%seg=();
	foreach $m (sort numerically keys %{$gnm{$j}} ) {				# m	= start
		$k4 = 0;	$line	= "";	$a	= 0;	$b	= "";
		foreach $n (sort numerically keys %{$gnm{$j}{$m}} ) {
			$gnm{$j}{$m}{$n}	=~ /\_(\d+)x*$/;	#	print STDERR $gnm{$j}{$m}{$n},"=name,cpN=$1\n";
			if ($k4 <= $1) {
				if ($k4 < $1) {
					$k4	= $1;						# copy number
					$a	= $n;						# end
					if ($b ne "") {$line	.= $b.",";}
					$b	= $gnm{$j}{$m}{$n};			# small RNA name
				} elsif ($k4 == $1 && $a < $n) {
					$a	= $n;
					if ($b ne "") {$line	.= $b.",";}
					$b	= $gnm{$j}{$m}{$n};
				}
			} else {
				$line	.= $gnm{$j}{$m}{$n}.",";
			}
	#		print STDERR "cpNum=$k4,end=$a,name=$b,all=$line\n";
		}
		$n	= $a;	$line	= $line.$b;	#	print STDERR "cpNum=$k4,end=$a,name=$b,all=$line\n\n";
#		foreach $n (sort numerically keys %{$gnm{$j}{$m}} ) {		# n	= end
			$a = $m > $endLen ? $m-$endLen : 0;
			$b = $a + $maxLen > $gnmLen[$k1]? $gnmLen[$k1]-$a : $maxLen;
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
	$foldFile		= "__miRiG_$start.$foldFileNum.fa";		# $foldFileNum
	$foldFileOut	= "__miRiG_$start.$foldFileNum.fa.out";
	open(FOLDTMPF, ">$foldFile") || die("Can not open file: $foldFile\n");
	foreach $key (sort keys %seg) {
		print FOLDTMPF ">$key\t",$seg{$key}->[0],"\n",$seg{$key}->[1],"\n";	$k4++;
	}
	close(FOLDTMPF);
	$foldFileNum++;
#	%seg=();
	print STDERR "\nThere are $k4 segments in $j -------------------------------------------------------\n\n";
#	next;

	print STDERR "\tperl $path/qsub_nonFasta1.1.pl -p \"perl $path/$ctFlt\" "."\n";
	system("perl $path/$qsub $qsub_set -p \"perl $path/$ctFlt\" -i \"-i\" -m \"$fold_ctFlt\" -q $foldFile -o \"2> $foldFileOut\" 2> $foldFile.err");	# --mfold=*,*,5 
	$a	= -1;
	while ($a == -1) {
		open(NUMTMPF, "$foldFile.err") || die("Can not open file: $foldFile.err\n");
		while (<NUMTMPF>) {
			if ($_=~/There\sare\s(\d+)\sjobs/) {
				$a	= $1; last;
			}
		}
		close(NUMTMPF);
	}
	print STDERR	"\tThere are $a jobs which were submitted!\n";
	if ($a	== 0) {
		%seg=();
		system("rm  $foldFile*");
		next;
	}
	$k2	= 0;	$b	= "";	$k4	= "";	my $k5	= 0;
	for (;$k2 < $a ;) {
		$k2 = 0;	$b	= "";	$k5	= 0;
		sleep(30);
		system("cat $foldFileOut*put > $foldFileOut");
		open(NUMTMPF, "$foldFileOut") || die("Can not open file: $foldFileOut\n");
		while (<NUMTMPF>) {
			if ($_=~/(\d+)\ssequences\sremained\safter\sstructure\scheck/) {
				$b	.= "\t\t" . $_;	#	print STDERR "k5=$k5+$1=";
				$k5+=$1;			#	print STDERR "$k5----------------\n";
			}
			$k4	.= $_;
	#		print STDERR "$_";
			if ($_ =~/Running\sfrom\s\[\d+\-\d+\-\d+.+\]\sto\s\[/) {
				$k2++;	$b	.= "\t" . $_;	#	print STDERR "\t$_";
			}
		}
		close(NUMTMPF);
	}
#	print STDERR "$k4\n;
	print STDERR "$b";	$k6	+= $k5;
	print STDERR "\n\tThere are $k5 sequences remained------------------------------------------------\n\n";
#	system("ls -l __$foldFile*ctFlt.ct > $k2");
	print STDERR "\tcat __$foldFile*ctFlt.ct/stem/txt > $foldFile.ctFlt.ct/stem/txt\n";
	
	if ($a	< 2) {
#		system("cat $foldFile*ctFlt.ct > $foldFile.ctFlt.ct");
	} else {
#		print STDERR "\ta=$a\n";
		system("cat __$foldFile*ctFlt.ct > $foldFile.ctFlt.ct");
	#	print STDERR "\tcat __$foldFile*ctFlt.stem > $foldFile.ctFlt.stem\n";
		system("cat __$foldFile*ctFlt.stem > $foldFile.ctFlt.stem");
	#	print STDERR "\tcat __$foldFile*ctFlt.txt > $foldFile.ctFlt.txt\n";
		system("cat __$foldFile*ctFlt.txt > $foldFile.ctFlt.txt");
	}
#	print STDERR "\tcat $foldFileOut* > $foldFileOut\n";
#	system("cat $foldFileOut*put > $foldFileOut");

	print STDERR "\tperl $path/searchseq1.4s.pl\n";
	system("perl $path/$searchSeq -i $foldFile -l $foldFile.ctFlt.txt > $foldFile.ctFlt.fa 2> $foldFile.err");

	print STDERR "\tcat $foldFile.ctFlt.fa/ct/stem/txt >> $infile.miRiG.fa/ct/stem/txt\n";
	system("cat $foldFile.ctFlt.fa >> $infile.miRiG.fa");
	system("cat $foldFile.ctFlt.ct >> $infile.miRiG.ct");
	system("cat $foldFile.ctFlt.stem >> $infile.miRiG.stem");
	system("cat $foldFile.ctFlt.txt >> $infile.miRiG.ctFlt.txt");

	system("du -sh > $foldFile.tmp");
	open(NUMTMPF, "$foldFile.tmp") || die("Can not open file: $foldFile\n");
	$a	= <NUMTMPF>;
	close(NUMTMPF);
#		print STDERR "\twc=$a";
	$a	=~ /^(\d+)/;
	if ($1 > 200) {
		print STDERR "\nWrong! This script occupy too many space: $a";
		system("rm __$foldFile* $foldFile*");
		%seg=();
		sub_end_program();
	}

	system("rm  __$foldFile* $foldFile*"); #$foldFileOut*put $foldFile.tmp $foldFile* 
#		system("perl ");
#	} else {
#		system("hybrid-ss-min -s DAT $foldFile > $foldFileOut");
#	}
#	$foldFileNum++;
	%seg=();
#	sub_end_program();
}

print STDERR "There are $k3,$k6 original,remained segments\n";
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
	print "-c	the parameter set for do_fold_ctFlt4miRNAInGnm.pl\n";
	print "			eg: -c \"\"; default: $fold_ctFlt\n";
	print "-h	display this lines\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
#	print "$0 -i blast_m8.blastn -t 90\n";
	print "$0 -i blast_m8.blastn -g genome.fa -l $len -p $path -o $output -q $qsub_set -c $fold_ctFlt\n";
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

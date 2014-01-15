#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.12.15
# Description:	read the blast m8 result and cut the segments from the genome files, then do unafold, then do ctStrucFilter
#**************************
# Version: 2.0	do ctFilter for each qsub results. use do_fold_ctFlt4miRNAInGnm.pl
# Version: 4.0	do unafold for $foldNum sequences each time.
# Version: 4.5	the infile may be blast or soap result. (New: identify soap result)
#**************************
# e-mail:highlei@gmail.com

my $version="4.5";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:g:d:l:p:o:q:c:e:t:k:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $gnmFiles	= (defined $opt_g) ? $opt_g : "";
my $smallRNA	= (defined $opt_d) ? $opt_d : "";
my $len			= (defined $opt_l) ? $opt_l : "20,100,20,300";
my $path		= (defined $opt_p) ? $opt_p : "./";
my $output		= (defined $opt_o) ? $opt_o : 0;		# 0: output the .fa and .ct files which 
my $qsub_set	= (defined $opt_q) ? $opt_q : "-b \"\" -s 32" ; # set -s and -r
my $fold_ctFlt	= (defined $opt_c) ? $opt_c : "-m 4 -n 2";

my $h2p_set		= (defined $opt_e) ? $opt_e : "-n,5,-r,0.75";
#my $minHitsNum	= (defined $opt_n) ? $opt_n : 5;
#my $minRatio	= (defined $opt_r) ? $opt_r : 0.75;

my $b2h_set		= (defined $opt_t) ? $opt_t : "0,-1,10,100,0,0,100";
#my $hsp_1st		= (defined $opt_s) ? $opt_s : 0;	# default:0: output all hsp; 1: output 1st hsp. 
#my $hit_1st		= (defined $opt_t) ? $opt_t : -1;	# default:-1: output all hits.
#my $evalue			= (defined $opt_e) ? $opt_e : 10;          #BlAST 搜索的域值
#my $identity		= (defined $opt_y) ? $opt_y : 100;
#my $score			= (defined $opt_c) ? $opt_c : 0;
#my $len			= (defined $opt_l) ? $opt_l : 0;	# 0:whole length match; -1: 1 base shorter than whole length
#my $overlap_len	= (defined $opt_p) ? $opt_p : 100;

my $pickup_set	= (defined $opt_k) ? $opt_k : "";

if ($opt_h || $infile eq ""){
	usage();
}

my $qsub		= "qsub4miRNAInGnm.pl -0 0";
my $ctFlt		= "do_fold_ctFlt4miRNAInGnm2.0.pl -0 0";#4miRNAInGnm
my $searchSeq	= "searchseq1.4s.pl";
my $searchCT	= "searchCTfile.pl";
my $pickup		= "pickupFromh2p4miRNAInGnm.pl -0 0";

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
	printf STDERR ("\n i  %55s : %-25s","input blast(m=8)/soap result",$infile);#%45s
	printf STDERR ("\n g  %55s : %-25s","input genome file",$gnmFiles);
	printf STDERR ("\n d  %55s : %-25s","input small RNA database (formatdb)",$smallRNA);
	printf STDERR ("\n l  %55s : %-25s","input length",$len);
	printf STDERR ("\n p  %55s : %-25s","input the path of perl scripts",$path);#%45s
	printf STDERR ("\n o  %55s : %-25s","output format",$output);
	printf STDERR ("\n q  %55s : %-25s","parameter set for qsub_noFasta.pl",$qsub_set);
	printf STDERR ("\n c  %55s : %-25s","parameter set for do_fold_ctFlt4miRNAInGnm.pl",$fold_ctFlt);
	printf STDERR ("\n e  %55s : %-25s","parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n t  %55s : %-25s","parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
	printf STDERR ("\n k  %55s : %-25s","parameter set for pickupFromh2p4miRNAInGnm.pl",$pickup_set);
#	if($zero==1) {printf STDERR ("\n z  %45s : %-25s","output coverage region?","1");}
#	elsif ($zero == 0) {printf STDERR ("\n z  %45s : %-25s","concise output","0");}
#	else {printf STDERR ("\n z  %45s : %-25s","output all genome","2");}
	printf STDERR ("\n x  %55s","exit the program!");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n------------------------------------------------------------\n\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input the blast(m=8)/soap result:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "g") {print STDERR "please input the genome file:\n"; $gnmFiles	= <STDIN>;$gnmFiles	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "d") {print STDERR "please input small RNA database,should be formatdb:\n"; $smallRNA	= <STDIN>;$smallRNA	=~s/[\s|\t|\r|\n]+$//g;}
#	elsif($yesorno eq "z") {$zero		= ($zero+1)%3;}
	elsif($yesorno eq "l") {print STDERR "please input len(endLen,minLen,addLen,maxLen):\n";$len	= <STDIN>;$len	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "p") {print STDERR "please input the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl):\n";$path	= <STDIN>;$path	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "o") {print STDERR "please input output format:\n";$output	= <STDIN>;$output	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "q") {print STDERR "please input the parameter set for qsub_noFasta.pl:\n";$qsub_set	= <STDIN>;$qsub_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "c") {print STDERR "please input the parameter set for do_fold_ctFlt4miRNAInGnm.pl:\n";$fold_ctFlt	= <STDIN>;$fold_ctFlt	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "e") {print STDERR "please input the parameter set for blast2histogram:\n";$b2h_set	= <STDIN>;$b2h_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "t") {print STDERR "please input the parameter set for histogram2predict4miRNAInGnm.pl:\n";$h2p_set	= <STDIN>;$h2p_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "k") {print STDERR "please input the parameter set for pickupFromh2p4miRNAInGnm.pl:\n";$pickup_set	= <STDIN>;$pickup_set	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {print STDERR ("============================================================\n");exit(0);}
}

	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %55s : %-25s","input blast(m=8)/soap result",$infile);#%45s
	printf STDERR ("\n g  %55s : %-25s","input genome file",$gnmFiles);
	printf STDERR ("\n d  %55s : %-25s","input small RNA database (formatdb)",$smallRNA);
	printf STDERR ("\n l  %55s : %-25s","input length",$len);
	printf STDERR ("\n p  %55s : %-25s","input the path of perl scripts",$path);#%45s
	printf STDERR ("\n o  %55s : %-25s","output format",$output);
	printf STDERR ("\n q  %55s : %-25s","parameter set for qsub_noFasta.pl",$qsub_set);
	printf STDERR ("\n c  %55s : %-25s","parameter set for do_fold_ctFlt4miRNAInGnm.pl",$fold_ctFlt);
	printf STDERR ("\n e  %55s : %-25s","parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n t  %55s : %-25s","parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
	printf STDERR ("\n k  %55s : %-25s\n","parameter set for pickupFromh2p4miRNAInGnm.pl",$pickup_set);

if ($len =~ /(\d+)\,(\d+)\,(\d+),(\d+)/) {
	$endLen	= $1;	$minLen	= $2;
	$addLen	= $3;	$maxLen	= $4;
	print STDERR "endLen=$endLen\tminLen=$minLen\taddLen=$addLen\tmaxLen=$maxLen\n\n";
} else {
	print STDERR "please input the correct len!\n";
	usage();
}
#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";

#------------------------------------------------v3.1 formatdb-----------------------------------------------
	$a	= 0;
	if (-e "$smallRNA.nhr") {} else {$a	= 1};
	if (-e "$smallRNA.nin") {} else {$a	= 1};
	if (-e "$smallRNA.nsq") {} else {$a	= 1};
	if ($a	== 1) {
		system("formatdb -p F -i $smallRNA");
	}
#------------------------------------------------v3.1 formatdb-----------------------------------------------

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
	} elsif ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+|-]+)\s+(\S+)\s+(\d+)\s+/) {
		$j	= $8;
		$m	= $9;
		$n	= $9+$6-1;
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

print	STDERR "\nLoad Alignment $infile OK\t$k1\n";
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
	print STDERR "Load Genome $gnmFiles OK\t",$i+1,"\n";
	############################################ chromosome legth ######################################################
	$m	= 1;
	for ($k2=0; $k2 <= $i;$k2++) {
		$gnmLen[$k2]	= length($genome[$k2]);
#		print STDERR "\tthe ",$m,"th sequence: length = $gnmLen[$k2] bp\n";
		$m++;
	}
#	print STDERR "\n";
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
	print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tdo $k4 sequences in $j\n";# -------------------------------------------------------\n\n";
	if ($k4 ==0) {
		%seg=();
		system("rm $foldFile*");
		next;
	}

	print STDERR "\tBegin submit jobs: $path/$ctFlt\n";
	$a	= -1;
	do {
		system("perl $path/$qsub $qsub_set -p \"perl $path/$ctFlt\" -i \"-i\" -m \"$fold_ctFlt -p $path -d $smallRNA -t $b2h_set -e $h2p_set\" -q $foldFile -o \"2> $foldFileOut\" 2> $foldFile.err");	# --mfold=*,*,5 
		$a	= -1;
		while ($a == -1) {
			open(NUMTMPF, "$foldFile.err") || die("Can not open file: $foldFile.err\n");
			while (<NUMTMPF>) {
				if ($_=~/There\sare\s(\d+)\sjobs/) {
					$a	= $1; last;
				} elsif ($_=~/failed/) {
					$a	= -2; last;
				}
			}
			close(NUMTMPF);
			if ($a == -2) {
				last;
			}
		}
		if ($a == -2) {
			print STDERR "Submit jobs failed!\n";
			print STDERR "Please re-run (no) or kill manually the successfully submited jobs and re-submit jobs (yes)\n";
			$k2=<STDIN>;
			if ($k2 =~/^n/) {
				system("rm  __$foldFile* $foldFile*");
				sub_end_program();
			}
		}
	} while ($a == -2) ;
	print STDERR	"\t$a jobs are submitted!\n";
	if ($a	== 0) {
		%seg=();
		system("rm __$foldFile* $foldFile*");
		next;
	}
	$k2	= 0;	$b	= "";	$k4	= "";	my $k5	= 0;	my $k7	= 0;	my $k8 = 0;
	for (;$k2 < $a ;) {
		$k2 = 0;	$b	= "";	$k4	= "";	$k5	= 0;	$k7	= 0;	$k8 = 0;
		sleep(30);
		system("cat $foldFileOut*put > $foldFileOut");
		open(NUMTMPF, "$foldFileOut") || die("Can not open file: $foldFileOut\n");
		while (<NUMTMPF>) {
			if ($_=~/(\d+)\ssequences\sremained\safter\sstructure\scheck/) {
	#			$b	.= "\t\t" . $_;	#	print STDERR "k5=$k5+$1=";
				$k5+=$1;			#	print STDERR "$k5----------------\n";
			}
			if ($_=~/(\d+)\ssequences which are deleted due to free energy > (\S+)/) {
	#			$b	.= "\t\t" . $_;	#	print STDERR "k5=$k5+$1=";
				$k7 += $1;			#	print STDERR "$k5----------------\n";
				$k8	= $2;
			}
			$k4	.= $_;
	#		print STDERR "$_";
			if ($_ =~/Running\sfrom\s\[\d+\-\d+\-\d+.+\]\sto\s\[/) {
				$k2++;
	#			$b	.= "\t" . $_;	#	print STDERR "\t$_";
			}
		}
		close(NUMTMPF);
	}
#	print STDERR "$k4\n;
#	print STDERR "$b";
	$k6	+= $k5;
	print STDERR "\tdelete $k7 sequences due to free energy > $k8\n";
	print STDERR "\tremain $k5 sequences after Structure check.\n";
#	system("ls -l __$foldFile*ctFlt.ct > $k2");
#	print STDERR "\tcat __$foldFile*ctFlt.ct/stem/txt/hist/h2p > $foldFile.ctFlt.ct/stem/txt/hist/h2p\n";
	
	if ($a	< 2) {
#		system("cat $foldFile*ctFlt.ct > $foldFile.ctFlt.ct");
	} else {
#		print STDERR "\ta=$a\n";
		system("cat __$foldFile*ctFlt.ct > $foldFile.ctFlt.ct");
	#	print STDERR "\tcat __$foldFile*ctFlt.stem > $foldFile.ctFlt.stem\n";
		system("cat __$foldFile*ctFlt.stem > $foldFile.ctFlt.stem");
	#	print STDERR "\tcat __$foldFile*ctFlt.txt > $foldFile.ctFlt.txt\n";
		system("cat __$foldFile*ctFlt.txt > $foldFile.ctFlt.txt");
		system("cat __$foldFile*ctFlt.hist > $foldFile.ctFlt.hist");
		system("cat __$foldFile*ctFlt.h2p > $foldFile.ctFlt.h2p");

	}
	$a	= `cat $foldFile.ctFlt.h2p | wc`;
#	$a	=~s/^[\s|\t]+//g;
	$a	=~ /^\s*(\d+)/;
	print STDERR "\tremain $1 after expression level.\n"; #$h2p_set
#	print STDERR "\tcat $foldFileOut* > $foldFileOut\n";
#	system("cat $foldFileOut*put > $foldFileOut");

#	print STDERR "\tperl $path/searchseq1.4s.pl\n";
	system("perl $path/$searchSeq -i $foldFile -l $foldFile.ctFlt.txt > $foldFile.ctFlt.fa 2> $foldFile.err");
	system("perl $path/$searchSeq -i $foldFile.ctFlt.hist -l $foldFile.ctFlt.txt > $foldFile.ctFlt.hist_ 2> $foldFile.err");
	system("perl $path/$searchSeq -i $foldFile.ctFlt.fa -l $foldFile.ctFlt.h2p > $foldFile.ctFlt.h2p.fa 2> $foldFile.err");

	system("perl $path/$pickup $pickup_set -i $foldFile.ctFlt.h2p > $foldFile.ctFlt.pkup 2> $foldFile.err");
	system("perl $path/$searchSeq -i $foldFile.ctFlt.fa -l $foldFile.ctFlt.pkup > $foldFile.ctFlt.pkup.fa 2> $foldFile.err");
	$a	= `awk \'/>/\' $foldFile.ctFlt.pkup.fa | wc`;
	$a	=~ /^\s*(\d+)/;
	print STDERR "\tremain $1 after removing overlap.\n"; 
#	print STDERR "\tcat $foldFile.ctFlt.fa/ct/stem/txt >> $infile.miRiG.fa/ct/stem/txt\n";
	system("cat $foldFile.ctFlt.fa >> $infile.miRiG.fa");
	system("cat $foldFile.ctFlt.ct >> $infile.miRiG.ct");
	system("cat $foldFileOut >> $infile.miRiG.err");
	system("cat $foldFile.ctFlt.stem >> $infile.miRiG.stem");
	system("cat $foldFile.ctFlt.txt >> $infile.miRiG.ctFlt.txt");
	system("cat $foldFile.ctFlt.hist_ >> $infile.miRiG.ctFlt.hist");
	system("cat $foldFile.ctFlt.h2p >> $infile.miRiG.ctFlt.h2p");
	system("cat $foldFile.ctFlt.h2p.fa >> $infile.miRiG.ctFlt.h2p.fa");
	system("cat $foldFile.ctFlt.pkup >> $infile.miRiG.ctFlt.pkup");
	system("cat $foldFile.ctFlt.pkup.fa >> $infile.miRiG.ctFlt.pkup.fa");

	$a	= `du -sh`;
#	print STDERR "\tspace=$a\n";
	$a	=~ /(\d+)G\s+\./;
#	print STDERR "\tspace=$1 G\n";
#	system("du -sh > $foldFile.tmp");
#	open(NUMTMPF, "$foldFile.tmp") || die("Can not open file: $foldFile\n");
#	$a	= <NUMTMPF>;
#	close(NUMTMPF);
#	$a	=~ /^(\d+)/;
#	print STDERR "\tspace=$1 G\n";
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
#	print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tfinish $j.\n";
	%seg=();
#	sub_end_program();
}

print STDERR "\nTotal sequences:\t$k3\n";
print STDERR "After Structure Check:\t$k6\n";
$a	= `cat $infile.miRiG.ctFlt.h2p | wc`;
$a	=~s/^[\s|\t]+//g;	$a	=~ /^(\d+)/;
print STDERR "After Expression Check:\t$1\n";
$a	= `awk \'/>/\' $infile.miRiG.ctFlt.pkup.fa | wc`;
$a	=~ /^\s*(\d+)/;
print STDERR "After Removing Overlap:\t$1\n"; 
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
	print "\t-i	<str>	input the  blast(m=8)/soap result file.";
	print " eg: blast_m8.blastn/soap.M0r2v0.soap\n";
	print "\t-g	<str>	input the genome file.";
	print " eg: rice.genome.fa\n";
	print "\t-d	<str>	input the small RNA database,should be formatdb.";
	print " eg: smallRNA.fa\n";
	print "\t-l	<int,int,int,int>	input len(endLen,minLen,addLen,maxLen).";
	print " [$len]\n";
	print "\t-p	<str>	the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl).";
	print " eg: ../soft; [$path] (current directory)\n";
	print "\t-o	<str>	output format.";
	print " [$output]\n";
	print "\t-q	<str>	the parameter set for qsub_noFasta.pl.";
	print " [$qsub_set]\n";
	print "\t-c	<str>	the parameter set for do_fold_ctFlt4miRNAInGnm.pl.";
	print " [$fold_ctFlt]\n";
	print "\t-e	<str>	the parameter set for blast2histogram.pl.";
	print " [$b2h_set]\n";
	print "\t-t	<str>	the parameter set for histogram2predict4miRNAInGnm.pl.";
	print " [$h2p_set]\n";
	print "\t-k	<str>	the parameter set for pickupFromh2p4miRNAInGnm.pl.";
	print " [$pickup_set]\n";
	print "\n\t-h	display this help\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
	print "$0 -i blast_m8.blastn/soap.M0r2v0.soap -g genome.fa -d smallRNA.fa -l $len -p $path -o $output -q \"-b \\\"\\\" -s 32\" -c \"$fold_ctFlt\"\n";
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
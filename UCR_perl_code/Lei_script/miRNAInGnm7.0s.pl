#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2010.09.21
# Description:	read the blast m8 result and cut the segments from the genome files, then do unafold, then do ctStrucFilter
#**************************
# Version: 2.0	do ctFilter for each qsub results. use do_fold_ctFlt4miRNAInGnm.pl
# Version: 4.0	do unafold for $foldNum sequences each time.
# Version: 4.5	the infile may be blast or soap result. (New: identify soap result)
# Version: 4.5s	swinlen	= length of smallRNA in do_fold_ctFlt4miRNAInGnm2.0s.pl
# Version: 5.0s	change the blast to soap in do_fold_ctFlt4miRNAInGnm3.0s.pl
# Version: 5.2s	redo the cycle if UNAfold is wrong. maxtime for redoing is 10.
# Version: 5.3s	do more than one chromosome each time.	Use $max_seq to control
# Version: 5.4s add parameter $start_from; it will be used to continue a uncompleted job. It assigns which chromosome this run start from. 
#				$N_num, delete the segments which contain more than $N_num Ns.
# Version: 5.5s use "undef(%hash);my %hash;" to replace "%hash=();"
# Version: 5.6s fix a bug in sub_readCTfile in do_fold_ctFlt4miRNAInGnm3.6s.pl
# Version: 6.0s use soap2p to filter h2p in do_fold_ctFlt4miRNAInGnm4.0s.pl
# Version: 6.1s fix soap2p, define $extent_len in do_fold_ctFlt4miRNAInGnm4.1s.pl
# Version: 6.2s fix a bug in cut genome
# Version: 6.3s delete the ct from %ct if fail in soap_check in do_fold_ctFlt4miRNAInGnm4.2s.pl.
# Version: 6.5s	plus/minus
# Version: 6.6s	use pickupFromh2p4miRNAInGnm6.6s.pl
# Version: 7.0s	do soap on genome. then use this soap result
#**************************
# e-mail:highlei@gmail.com

my $version="7.0s";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:g:d:l:p:o:s:q:c:a:j:m:f:e:t:k:S:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $gnmFiles	= (defined $opt_g) ? $opt_g : "";
my $smallRNA	= (defined $opt_d) ? $opt_d : "";
my $len			= (defined $opt_l) ? $opt_l : "20,100,20,300";
my $path		= (defined $opt_p) ? $opt_p : "./";
my $output		= (defined $opt_o) ? $opt_o : 0;		# 0: output the .fa and .ct files which 
my $splitQuery	= (defined $opt_s) ? $opt_s : 4;
my $qsub_set	= (defined $opt_q) ? $opt_q : "-b \"\"" ; # set -s and -r
my $fold_ctFlt	= (defined $opt_c) ? $opt_c : "-m 4 -n 1";
my $mismatch_in_gnm	= (defined $opt_a) ? $opt_a : 0;
my $blastORsoap	= (defined $opt_j) ? $opt_j : "blast";
my $max_seq		= (defined $opt_m) ? $opt_m : 10000;
my $start_from	= (defined $opt_f) ? $opt_f - 1 : 0;

my $N_num		= (defined $opt_N) ? $opt_N : 25;

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
my $soap_result	= (defined $opt_S) ? $opt_S : "";



if ($opt_h || $infile eq ""){
	usage();
}

my $qsub		= "qsub4miRNAInGnm1.1.pl -0 0";
my $ctFlt		= "do_fold_ctFlt4miRNAInGnm$version.pl -0 0";#4miRNAInGnm
my $searchSeq	= "searchseq1.7s.pl";
my $searchCT	= "searchCTfile.pl";
my $pickup		= "pickupFromh2p4miRNAInGnm$version.pl -0 0";
my $searchLine	= "searchLineACList1.2s.pl";
sub numerically{$a<=>$b};

use FileHandle;
use strict;



my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$k4,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@genome,@gnmName,@gnmLen);
my (%gnm,%seg,%strand);
my $key="";
my ($endLen,$minLen,$addLen,$maxLen,$foldFileNum,$soap_M);
my ($foldFile,$foldFileOut,$soapOut);

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
	printf STDERR ("\n j  %55s : %-25s","use blast or soap",$blastORsoap);
	printf STDERR ("\n s  %55s : %-25s","split query",$splitQuery);
	printf STDERR ("\n m  %55s : %-25s","max sequence number for one submission",$max_seq);
	printf STDERR ("\n f  %55s : %-25s","from which chromosome this run starts",$start_from);
	printf STDERR ("\n a  %55s : %-25s","soap M (mismatch num when map)",$mismatch_in_gnm);
	printf STDERR ("\n S  %55s : %-25s","soap result file (-d VS -g)",$soap_result);
	printf STDERR ("\n q  %55s : %-25s","parameter set for qsub_noFasta.pl",$qsub_set);
	printf STDERR ("\n c  %55s : %-25s","parameter set for do_fold_ctFlt4miRNAInGnm.pl",$fold_ctFlt);
	printf STDERR ("\n t  %55s : %-25s","parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n e  %55s : %-25s","parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
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
	elsif($yesorno eq "j") {print STDERR "please input use blast or soap:\n";$blastORsoap	= <STDIN>;$blastORsoap	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input how many files split query to:\n";$splitQuery	= <STDIN>;$splitQuery	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "m") {print STDERR "please input max sequence number for one submission:\n";$max_seq	= <STDIN>;$max_seq	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input from which chromosome this run starts:\n";$start_from	= <STDIN>;$start_from	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "a") {print STDERR "please input soap M (mismatch num when map):\n";$mismatch_in_gnm	= <STDIN>;$mismatch_in_gnm	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "S") {print STDERR "please input soap result file (-d VS -g):\n";$soap_result	= <STDIN>;$soap_result	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "q") {print STDERR "please input the parameter set for qsub_noFasta.pl:\n";$qsub_set	= <STDIN>;$qsub_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "c") {print STDERR "please input the parameter set for do_fold_ctFlt4miRNAInGnm.pl:\n";$fold_ctFlt	= <STDIN>;$fold_ctFlt	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "t") {print STDERR "please input the parameter set for blast2histogram:\n";$b2h_set	= <STDIN>;$b2h_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "e") {print STDERR "please input the parameter set for histogram2predict4miRNAInGnm.pl:\n";$h2p_set	= <STDIN>;$h2p_set	=~s/[\s|\t|\r|\n]+$//g;}
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
	printf STDERR ("\n s  %55s : %-25s","split query",$splitQuery);
	printf STDERR ("\n m  %55s : %-25s","max sequence number for one submission",$max_seq);
	printf STDERR ("\n f  %55s : %-25s","from which chromosome this run starts",$start_from);
	printf STDERR ("\n a  %55s : %-25s","soap M (mismatch num when map)",$mismatch_in_gnm);
	printf STDERR ("\n S  %55s : %-25s","soap result file (-d VS -g)",$soap_result);
	printf STDERR ("\n q  %55s : %-25s","parameter set for qsub_noFasta.pl",$qsub_set);
	printf STDERR ("\n c  %55s : %-25s","parameter set for do_fold_ctFlt4miRNAInGnm.pl",$fold_ctFlt);
	printf STDERR ("\n t  %55s : %-25s","parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n e  %55s : %-25s","parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
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

if ($mismatch_in_gnm == 4 || $mismatch_in_gnm == 2) {
	$soap_M	= 4;
	$mismatch_in_gnm	= 2;
} else {
	$soap_M	= $mismatch_in_gnm;
}
#------------------------------------------------v3.1 formatdb-----------------------------------------------
	$a	= 0;
	if (-e "$smallRNA.nhr") {} else {$a	= 1};
	if (-e "$smallRNA.nin") {} else {$a	= 1};
	if (-e "$smallRNA.nsq") {} else {$a	= 1};
	if ($a	== 1 && $blastORsoap =~/^b/) {
		system("formatdb -p F -i $smallRNA");
	}
#------------------------------------------------v3.1 formatdb-----------------------------------------------

if ($soap_result eq "") {
	system("2bwt-builder $gnmFiles");
	$soap_result	= "__miRiG_$start.sRNA.VS.gnm.soap";
	system("soap -D $gnmFiles.index -a $smallRNA -M $soap_M -r 2 -v 0 -o $soap_result 2> $soap_result.err");
}
############################################ read file ######################################################
$soapOut	= "__miRiG_$start.sRNA.VS.chr.soap";

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
	} elsif ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+|-]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+/) {
		if ($10 > $mismatch_in_gnm) {
			next;
		}
		$j	= $8;	$a	= $7;
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
				$strand{$i}{"$j.$m.$n"}	= $a;
			}													#	print STDERR $gnm{$j}{$m}{$n},"\n";
#			print STDERR "This occurs more than one time:",$_;
		} else {
			$gnm{$j}{$m}{$n}	= $1;
			$strand{$1}{"$j.$m.$n"}	= $a;
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
	print STDERR "Load Genome $gnmFiles OK\t",$i+1,"\n\n";
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

my $ESTorGNM= $i < 100 ? 0 : 1;
print STDERR "Genome number < 100\t",$ESTorGNM,"\n\n";

if ($start_from <= 0) {
	$k1 = `rm $infile.miRiG* 2>&1`;
	$k1 = 0;
} else {
	$k1	= $start_from;
}
$foldFileNum	= $k1;

$k3	= 0;	my $k6	= 0;	my $k10	= 0;	my $k11	= 0;	my	$k12	= 0;	my $cpn1= 0;	my $chrName="";
for (;$k1 <= $i ;$k1++) {
	$j	= $gnmName[$k1];	$k12++;
	if ($k11 > $max_seq) {
		%seg=();	$k11	= 0;	$k12	= 1;
	}
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
					$seg{$key}->[0]	=~ /\_(\d+)x*\;\d+\,\d+\;\d+\,\d+$/;
					$cpn1	= $1;
					$line	=~ /\_(\d+)x*$/;
					if ($cpn1 > $1) {
						$seg{$key}->[0]	= $line.";$m,$n;".($m-$a).",".($n-$a)."\t".$seg{$key}->[0];
					} else {
						$seg{$key}->[0]	.= "\t".$line.";$m,$n;".($m-$a).",".($n-$a);
					}
				} else {
#					$line	=~ /(\w+)$/;
#					if ($strand{$1}{"$j.$m.$n"} =~/\+/) {
						$seg{$key}->[0]	= $line.";$m,$n;".($m-$a).",".($n-$a);
						$seg{$key}->[1]	= substr($genome[$k1],$a,$k2);
#					} elsif ($strand{$1}{"$j.$m.$n"} =~/\-/) {
#						$seg{$key}->[0]	= $line.";$m,$n;".($k2-$n+$a+1).",".($k2-$m+$a+1);
#						$seg{$key}->[1]	= reverseDNAString(substr($genome[$k1],$a,$k2));
#					} else {
#						print STDERR $strand{$1}{"$j.$m.$n"},",$j,$m,$n\t$1\n";
#						die("Wrong in strand 333\n");
#					}
					$k3++;	$k11++;
				}
	#			print ">$j\_$a\_$k2\t",$gnm{$j}{$m}{$n},"\n",substr($genome[$k1],$a,$k2),"\n";
			}
			$a	= $n + $endLen > $gnmLen[$k1]? $gnmLen[$k1] : $n+$endLen;
			$b	= $a - $maxLen > 0 ? $maxLen : $a;
			for ($k2 = $minLen;$k2 <= $b ;$k2+=$addLen) {
				$key	= "$j\_" . ($a-$k2) . "\_$k2";
				if (exists($seg{$key})) {
					$seg{$key}->[0]	=~ /\_(\d+)x*\;\d+\,\d+\;\d+\,\d+$/;
					$cpn1	= $1;
					$line	=~ /\_(\d+)x*$/;
					if ($cpn1 > $1) {
						$seg{$key}->[0]	= $line.";$m,$n;".($m-$a+$k2).",".($n-$a+$k2)."\t".$seg{$key}->[0];
					} else {
						$seg{$key}->[0]	.= "\t".$line.";$m,$n;".($m-$a+$k2).",".($n-$a+$k2);
					}
				} else {
					$seg{$key}->[0]	= $line.";$m,$n;".($m-$a+$k2).",".($n-$a+$k2);
					$seg{$key}->[1]	= substr($genome[$k1],$a-$k2,$k2);	$k3++;	$k11++;
				}
	#			print ">$j\_",$a-$k2,"\_$k2\t",$gnm{$j}{$m}{$n},"\n",substr($genome[$k1],$a-$k2,$k2),"\n";
			}
#		}
	}
	foreach $key ( keys %seg) {
		if ($key !~ /^$j/) {
			next;
		}
		$seg{$key}->[0]	=~ /(\w+)\;(\d+)\,(\d+)\;\d+\,\d+$/;
#		print STDERR "key=$key,",$seg{$key}->[0],"\n";
		$line	= $1;	$m	= $2;	$n	= $3;
		if ($strand{$line}{"$j.$m.$n"} =~ /\+/) {
		} elsif ($strand{$line}{"$j.$m.$n"} =~ /\-/) {
#			print STDERR $seg{$key}->[1],"\n";
			$seg{$key}->[1]	= reverseDNAString($seg{$key}->[1]);
#			print STDERR $seg{$key}->[1],"\n";
#			$k2	= length($seg{$key}->[1]);
#			print STDERR "$key,",$seg{$key}->[0],"\n";
			$seg{$key}->[0]	=~ /(\d+)\,(\d+)$/;
#			$a	= $k2-$2+1;	$b	= $k2-$1+1;
			$seg{$key}->[0]	=~ s/(\d+)\,(\d+)$/-$1\,-$2/;#s/\d+\,\d+$/$a\,$b/;
#			print STDERR $seg{$key}->[0],"\n";
#			$k2=<STDIN>;
		} else {
			print STDERR $strand{$line}{"$j.$m.$n"},"\t",$seg{$key}->[0],"\t",$seg{$key}->[1],"\tkey=$key,$j,$m,$n\t$line\n";
			die("Wrong in strand 333\n");
		}
	}
#	undef(%strand);	my %strand=();
	$k4	= 0;
	$foldFile		= "__miRiG_$start.$foldFileNum.fa";		# $foldFileNum
	$foldFileOut	= "__miRiG_$start.$foldFileNum.fa.out";
	print STDERR "\t$k11 sequences till $j\n";
#	system("less $soap_result | gawk '{if(\$8==\"$j\"){print \$_}}' >> $soapOut");
#	$chrName	.= "$j\n";
	if($ESTorGNM==0){system("perl -e 'print \"$j\n\"' >> $foldFile.lis");}
#	system("perl -e 'print \"$j\n\"'");
	if ($k11	< $max_seq && $k1 < $i) {
		next;
	}
	open(FOLDTMPF, ">$foldFile") || die("Can not open file: $foldFile\n");
	foreach $key (sort keys %seg) {
		if ($seg{$key}->[1] =~ /N{$N_num}/i) {
			next;
		}
		print FOLDTMPF ">$key\t",$seg{$key}->[0],"\n",$seg{$key}->[1],"\n";	$k4++;
	}
	close(FOLDTMPF);
	$foldFileNum++;
	undef(%seg);	my %seg=();
	print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tdo $k4 sequences in $j\n";# -------------------------------------------------------\n\n";
	if ($k4 ==0) { # || $j eq "Chr1") {
		%seg=();	$chrName	= "";
		system("rm $foldFile* $soapOut");
		next;
	}
#	print STDERR "chr=$chrName;\t";
#	system("perl -e '{print \"$chrName\n\";}' >> $foldFile.lis");
	if ($ESTorGNM	== 0) {
		system("perl $path/$searchLine -i $soap_result -l $foldFile.lis -d 8 > $soapOut 2> $foldFile.lis.err");
	} else {
		system("cat $soap_result > $soapOut");
	}
	print STDERR "\tBegin submit jobs: $path/$ctFlt\n";	#$a=<STDIN>;
	$a	= -1;
	do {
		system("perl $path/$qsub $qsub_set -s $splitQuery -p \"perl $path/$ctFlt\" -i \"-i\" -m \"$fold_ctFlt -p $path -j $blastORsoap -S $soapOut -d $smallRNA -t $b2h_set -e $h2p_set\" -q $foldFile -o \"2> $foldFileOut\" 2> $foldFile.err");	# --mfold=*,*,5 
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
	$k2	= 0;	$b	= "";	$k4	= "";	my $k5	= 0;	my $k7	= 0;	my $k8 = 0;	my $k9	= 0;
#	my $redoFlag	= 0;
	for (;$k2 < $a ;) {
		$k2 = 0;	$b	= "";	$k4	= "";	$k5	= 0;	$k7	= 0;	$k8 = 0;
		sleep(30);
#		$redoFlag	= `cat $foldFileOut*put > $foldFileOut`;
		system("cat $foldFileOut*put > $foldFileOut 2>&1");
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
			if ($_ =~/Cannot\sfinish\sthis\sfold/) {
				$k9	= 1;	$k2++;
			}
			if ($_ =~/cat\:\s+$foldFileOut.+No such file or directory/) {
				$k9	= 1;	$k2++;
			}
		}
		close(NUMTMPF);
	}
	if ($k9	== 1) {
		$k1	-= $k12;	$k10++;
		system("rm __$foldFile* $foldFile*");
		%seg=();	$chrName	= "";
		if ($k10	> 10) {
			print STDERR "Error: cannot finish the cycle for $j!\n";
			sub_end_program();
		}
		print STDERR "[***| Redo the cycle ($k10) in $j |***]\n";
		next;
	}
	$k10	= 0;
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
		system("cat __$foldFile*ctFlt.spck > $foldFile.ctFlt.spck");
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

	system("perl $path/$pickup $pickup_set -i $foldFile.ctFlt.h2p -f $foldFile.ctFlt.h2p.fa > $foldFile.ctFlt.pkup 2> $foldFile.err");
	system("perl $path/$searchSeq -i $foldFile.ctFlt.fa -l $foldFile.ctFlt.pkup > $foldFile.ctFlt.pkup.fa 2> $foldFile.err");
	$a	= `awk \'/>/\' $foldFile.ctFlt.pkup.fa | wc`;
	$a	=~ /^\s*(\d+)/;
	print STDERR "\tremain $1 after removing overlap.\n\n"; 
#	print STDERR "\tcat $foldFile.ctFlt.fa/ct/stem/txt >> $infile.miRiG.fa/ct/stem/txt\n";
	system("cat $foldFile.ctFlt.fa >> $infile.miRiG.fa");
	system("cat $foldFile.ctFlt.ct >> $infile.miRiG.ct");
	system("cat $foldFileOut >> $infile.miRiG.err");
	system("cat $foldFile.ctFlt.stem >> $infile.miRiG.stem");
	system("cat $foldFile.ctFlt.txt >> $infile.miRiG.ctFlt.txt");
	system("cat $foldFile.ctFlt.spck >> $infile.miRiG.ctFlt.spck");
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
		%seg=();	$chrName	= "";
		sub_end_program();
	}
	system("rm  __$foldFile* $foldFile* $soapOut"); #$foldFileOut*put $foldFile.tmp $foldFile* 
#		system("perl ");
#	} else {
#		system("hybrid-ss-min -s DAT $foldFile > $foldFileOut");
#	}
#	$foldFileNum++;
#	print STDERR "\nNow = [",sub_format_datetime(localtime(time())),"]\tfinish $j.\n";
	%seg=();	$chrName	= "";
#	sub_end_program();
}
if ($soap_result	eq "__miRiG_$start.sRNA.VS.gnm.soap") {
	system("rm $soap_result*");
}
print STDERR "\nTotal sequences:\t$k3\n";
print STDERR "After Structure Check:\t$k6\n";
$a	= `cat $infile.miRiG.ctFlt.h2p | wc`;
$a	=~s/^[\s|\t]+//g;	$a	=~ /^(\d+)/;
print STDERR "After Expression Check:\t$1\n";
$a	= `perl -e '\$i=0;while(<>){\$i++ if /^>/} print \"\$i\n\"' $infile.miRiG.ctFlt.pkup.fa`;
print STDERR "After Removing Overlap:\t$a\n"; 
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
	print "\t-p	<str>	the path of perl scripts.";
	print " eg: ../soft; [$path] (current directory)\n";
#	print "\t\t\tqsub4miRNAInGnm.pl,\n\t\t\tdo_fold_ctFlt4miRNAInGnm3.0s.pl,";
#	print "\n\t\t\tsearchseq1.4s.pl,\n\t\t\tsearchCTfile.pl,\n\t\t\tblast2histogram.pl,";
#	print "\n\t\t\thistogram2predict4miRNAInGnm.pl,\n\t\t\tpickupFromh2p4miRNAInGnm.pl\n";
	print "\t-j	<str>	use blast or soap.";
	print " [$blastORsoap]\n";
	print "\t-s	<int>	the number to split the query into.";
	print " [$splitQuery]\n";
	print "\t-m	<int>	max sequence number for one submission.";
	print " [$max_seq]\n";
	print "\t-f	<int>	from which chromosome this run starts.";
	print " [$start_from]\n";
	print "\t-a	<int>	soap M (mismatch num when map).";
	print " [$mismatch_in_gnm]\n";
	print "\t-S	<int>	soap result file (-d VS -g).";
	print " [$soap_result]\n";
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
	print "$0 -i blast_m8.blastn/soap.M0r2v0.soap -g genome.fa -d smallRNA.fa -l $len -j $blastORsoap -p $path -q \"$qsub_set\" -s $splitQuery -c \"$fold_ctFlt\"\n";
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
############################################################################################################
######################                  reverseDNAString
############################################################################################################

sub reverseDNAString 
{
	my($rdstr)	= @_;
	my ($sr1,$sr2);
	$rdstr	= reverse($rdstr);
	$rdstr	=~tr/ACGTRYMKacgtrymk/TGCAYRKMtgcayrkm/;
	return $rdstr;
}
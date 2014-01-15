#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			histogram2predict4miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.03
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2010.09.21
# Description:	read blast2histogram.pl result and ctStructure2.1.pl result and predict the possible miRNAs.
#**************************
# Version: 2.0		deal with the multi-structure of one sequence (hybrid-ss-min --mfold=*,*,5 will output 5 max optimal structures for one sequence)
# Version: 3.0		do each chromosome
# Version: 6.5s		plus/minus
# Version: 6.6s		use pickupFromh2p4miRNAInGnm6.6s.pl
# Version: 7.0s		do soap on genome. then use this soap result
#**************************
#refer to similarity4phy1.4.pl

# e-mail:highlei@gmail.com

my $version="7.0s";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:c:n:r:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $histFile	= $opt_i;
my $ctFltFile	= $opt_c;
my $minHitsNum	= (defined $opt_n) ? $opt_n : 5;
my $minRatio	= (defined $opt_r) ? $opt_r : 0.75;
my $winLen		= (defined $opt_w) ? $opt_w : 21;


if ($opt_h || $histFile eq "" || $ctFltFile eq "") { #|| $batchfile eq "" || $output eq "") {
	usage();
}
sub numerically{$a<=>$b};
#sub sub_slideWindow;
#sub sub_numOfDiff;
my (%ctFlt,%mrna,%chr);

use FileHandle;
use strict;



my ($i,$j,$k,$num,$len,$ntlen,$k1,$k2,$k3,$k4,$m,$n,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@buf,@bufo,@genome,@gnmName,);

my $key="";
my $key2	="";
my $flag=0;
my ($opNum,$fore,$back,$gnmNum,$site);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %40s : %-25s","input blast2histogram.pl result file name",$histFile);#%45s
	printf STDERR ("\n c  %40s : %-25s","input ctStructure2.1.pl result file name",$ctFltFile);#%45s
	printf STDERR ("\n n  %40s : %-25s","input minnimum hits number",$minHitsNum);
	printf STDERR ("\n r  %40s : %-25s","input minnimum ratio of duplex",$minRatio);
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input blast2histogram.pl result file name:\n"; $histFile	= <STDIN>;	$histFile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "c") {print STDERR "please input ctStructure2.1.pl result file name:\n"; $ctFltFile	= <STDIN>;	$ctFltFile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "n") {print STDERR "please input minnimum hits number\n"; $minHitsNum	= <STDIN>;$minHitsNum	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "r") {print STDERR "please input minimum ratio of duplex in segments:\n";$minRatio	= <STDIN>;$minRatio	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "x") {sub_end_program();;exit(0);}
}


#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";
############################################ read histFile files ######################################################
$i	= -1;
$file = new FileHandle ("$histFile") || die("Cannot open file:$histFile\n");
while (<$file>) {
	$_=~s/^[\s|\t]+//g;
	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^>(\S+)/) {
	#	if (exists($genome{$1})) {
	#		die("There are something wrong in $histFile:$1\n");
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
print STDERR "There are ",$i+1," sequences in file $histFile\n";
$gnmNum	= $i+1;
############################################ read ctFltFile files ######################################################
$i	= -1;	$n	= 0;
$file = new FileHandle ("$ctFltFile") || die("Cannot open file:$ctFltFile\n");
while (<$file>) {
	if ($_ =~/^>?([^\_]+)\_/) {
		if (exists($chr{$1}{$_})) {
			print STDERR "multi times:$_";
		} else {
			$chr{$1}{$_}	= 1;
		}
		$i++;
	} else {
		print STDERR "match wrong in 106:$_\n";
	}
}
close $file || die("Wrong!");
print STDERR "There are ",$i+1," sequences in file $ctFltFile\n";
############################################ cut genome and do fold ######################################################

$i	= -1;	$n	= 0;	my $ch;
my ($a1,$a2,$a3,$a4);
foreach $ch (sort keys %chr) {
	%ctFlt	= ();	print STDERR "$ch\n";
	foreach $line (sort keys %{$chr{$ch}}) {	#	print STDERR $line;
		if ($line =~/^>?(\S+)\s+(.+)$/) {
			$key	= $1;	$j	= $2;	splice(@buf,0);		#	print STDERR $key," !\n";
			if (!exists($ctFlt{$key})) {
				$ctFlt{$key}{"num"}		= 1;	$n	= 1;
			} else {
				$ctFlt{$key}{"num"}++;	$n	= $ctFlt{$key}{"num"};	#print STDERR "n=$n,num=",$ctFlt{$key}{"num"},"\n";
			}
			$key	=~/_(\d+)$/;
			my $mLen	= $1;
			@buf    = split(/\;\s?/,$j);
			for ($k1 = 0;$k1 < @buf ;$k1++) {
				if ($buf[$k1] =~/(\d+)\.\.(\d+)\,(\d+)\.\.(\d+)/) {
					for ($k2 = $1;$k2 <= $2 ;$k2++) {
						$ctFlt{$key}{$n}{$buf[$k1]}{"L"}{$k2}	= 0;
					}
					if ($3 > $4) {
						for ($k2 = $3;$k2 >= $4 ;$k2--) {
							$ctFlt{$key}{$n}{$buf[$k1]}{"R"}{$k2}	= 0;
						}
					} else {
						for ($k2 = $3;$k2 <= $4 ;$k2++) {
							$ctFlt{$key}{$n}{$buf[$k1]}{"R"}{$k2}	= 0;
						}
					}
				} elsif ($buf[$k1] =~/(-\d+)\.\.(-\d+)\,(-\d+)\.\.(-\d+)/) {
					$a1	= $2+$mLen;	$a2 = $1+$mLen;	$a3 = $4+$mLen;	$a4 = $3+$mLen;
					for ($k2 = $a1;$k2 <= $a2 ;$k2++) {
						$ctFlt{$key}{$n}{$buf[$k1]}{"L"}{$k2}	= 0;
					}
					if ($a3 > $a4) {
						for ($k2 = $a3;$k2 >= $a4 ;$k2--) {
							$ctFlt{$key}{$n}{$buf[$k1]}{"R"}{$k2}	= 0;
						}
					} else {
						for ($k2 = $a3;$k2 <= $a4 ;$k2++) {
							$ctFlt{$key}{$n}{$buf[$k1]}{"R"}{$k2}	= 0;
						}
					}
				} else {
					print STDERR "match wrong in 101:$buf[$k1]\n";
				}
			}
			$i++;
		} else {
			print STDERR "match wrong in 106:$_\n";
		}
	}
	print STDERR "\tThere are $i lines in $ch\n";
#----------------------------------------------------------------------------------------------------------------

for ($k1 = 0;$k1 <= $gnmNum ;$k1++) {	
#	print STDERR "k1=$k1, $gnmName[$k1] begin!\n";
	if (!exists($ctFlt{$gnmName[$k1]})) {
#		print STDERR $gnmName[$k1]," no ct structure info!\n";
		next;
	}			#	print STDERR $gnmName[$k1]," start!\n";
	splice(@buf,0);	$k	= 0;							# k the sum of all histogram.
	@buf    = split(/\s+/,$genome[$k1]);
	for ($k2 = 0;$k2 < @buf ;$k2++) {
		$k	+= $buf[$k2];
	}			#	print STDERR $gnmName[$k1]," cycle!\n";
	foreach $key (sort keys %{$ctFlt{$gnmName[$k1]}}) {		#print STDERR $gnmName[$k1]," in cycle!\n";
		if ($key eq "num") {
			next;
		}
		%mrna	= ();	$k4 = 0;		#	print STDERR $gnmName[$k1]," in mrna!\n";
		foreach $key2 (sort keys %{$ctFlt{$gnmName[$k1]}{$key}}) {
			$m	= 0;	$a	= 0;	$k3	= 0;	$b	= 0;	$n	= 0;
			foreach $k2 (sort numerically keys %{$ctFlt{$gnmName[$k1]}{$key}{$key2}{"L"}}) {
				if ($b	== 0) {
					$k3	= $k2;
					$b++;
				}
				$ctFlt{$gnmName[$k1]}{$key}{$key2}{"L"}{$k2}	= $buf[$k2-1];
				$m	+= $buf[$k2-1];	$n	= $k2;
				$a	= $a > $buf[$k2-1] ? $a : $buf[$k2-1];
			}
#			print STDERR "m=$m,L=$k3..$n\t";
			$mrna{$k4}{"L"}{"start"}	= $k3;
			$mrna{$k4}{"L"}{"end"}		= $n;
			$mrna{$k4}{"L"}{"num"}		= $m;
			$mrna{$k4}{"L"}{"hit"}		= $a;
			$mrna{$k4}{"L"}{"site"}		= $key2;
		#	$ctFlt{$gnmName[$k1]}{$key}{"L-ratio"}	= 1.0*$m/$k;
			$m	= 0;	$a	= 0;	$k3	= 0;	$b	= 0;
			foreach $k2 (sort numerically keys %{$ctFlt{$gnmName[$k1]}{$key}{$key2}{"R"}}) {
				if ($b	== 0) {
					$k3	= $k2;
					$b++;
				}
				$ctFlt{$gnmName[$k1]}{$key}{$key2}{"R"}{$k2}	= $buf[$k2-1];
				$m	+= $buf[$k2-1];	$n	= $k2;
				$a	= $a > $buf[$k2-1] ? $a : $buf[$k2-1];
			}
#			print STDERR "m=$m,R=$k3..$n\n";
			if ($mrna{$k4}{"L"}{"start"} > $k3 && $mrna{$k4}{"L"}{"end"} > $n) {
				$mrna{$k4}{"R"}{"start"}	= $mrna{$k4}{"L"}{"start"};
				$mrna{$k4}{"R"}{"end"}		= $mrna{$k4}{"L"}{"end"};
				$mrna{$k4}{"R"}{"num"}		= $mrna{$k4}{"L"}{"num"};
				$mrna{$k4}{"R"}{"hit"}		= $mrna{$k4}{"L"}{"hit"};
				$mrna{$k4}{"L"}{"start"}	= $k3;
				$mrna{$k4}{"L"}{"end"}		= $n;
				$mrna{$k4}{"L"}{"num"}		= $m;
				$mrna{$k4}{"L"}{"hit"}		= $a;
			} else {
				$mrna{$k4}{"R"}{"start"}	= $k3;
				$mrna{$k4}{"R"}{"end"}		= $n;
				$mrna{$k4}{"R"}{"num"}		= $m;
				$mrna{$k4}{"R"}{"hit"}		= $a;
			}
#			print STDERR "Lsite: sn=",$mrna{$k4}{"L"}{"site"},"\tnum:L=",$mrna{$k4}{"L"}{"num"},"\tR=",$mrna{$k4}{"R"}{"num"},"\n";
			$k4++;
		}							#	print STDERR $gnmName[$k1]," before sub_link_max!\n";
		do {
			$m	= sub_link_max();
		} while ($m > 0);
		$m	= 0;	$n	= 0;	$a	= 0;	$b	= "";
		foreach $k4 (sort numerically keys %mrna) {
			if ($mrna{$k4}{"L"}{"hit"} < $minHitsNum && $mrna{$k4}{"R"}{"hit"} < $minHitsNum) {
				next;
			}
#			print STDERR "$gnmName[$k1],",$mrna{$k4}{"L"}{"start"},"..",$mrna{$k4}{"L"}{"end"},";",$mrna{$k4}{"R"}{"start"},"..",$mrna{$k4}{"R"}{"end"},"\n";
	##		$m	+= $mrna{$k4}{"L"}{"hit"}*$winLen;	$n	+= $mrna{$k4}{"R"}{"hit"}*$winLen;
			$m	+= $mrna{$k4}{"L"}{"num"};	$n	+= $mrna{$k4}{"R"}{"num"};
			$a++;	$b	.= $mrna{$k4}{"L"}{"site"}.";";
		}							#	print STDERR $gnmName[$k1]," ratio\n";
		if ($a == 0 || 1.0*($m+$n)/$k < $minRatio) {
			next;
		}
		print $gnmName[$k1],"\t",$a,";$b\t";	#	print STDERR $gnmName[$k1]," output!\n";
		printf("%.3f\t%.3f\t%.3f\n",1.0*$m/$k,1.0*$n/$k,1.0*($m+$n)/$k);
	}
}
#-------------------------------------------------------------------------------------------------------

}



############################################ output ##########################################################

sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################


sub usage
{
	print "Contact : Gaolei <highlei\@highlei.com>";
	print "\nProgram : $0\nVersion: $version\n";
	print "retro-translate amino-acid sequences (alignment) to nucleotide sequences\n";
	print "Usage:\n	$0 \n";
	print "-i	input blast2histogram.pl result file name\n";
	print "			eg: blastHist.txt\n";
	print "-c	input ctStructure2.1.pl result file name\n";
	print "			eg: ctFilter2.1.txt\n";
	print "-n	input minnimum hits number\n";
	print "			eg: -n 30; default: $minHitsNum\n";
	print "-r	input minimum ratio of duplex in segments\n";
	print "			eg: -r 0.5; default: $minRatio\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i blastHist.txt -c ctFilter2.1.txt -n 10 -r 0.75\n";# -m \"-k 3 -s DAT\"\n";
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

############################################################################################################
######################                  sub_link_max
############################################################################################################

sub sub_link_max
{
	my ($s00,$s01,$s10,$s11,$s1,$sk,$sf,$sf2);
	my ($schr,$sn,$sn1,$sleft,$sflag,$sflag2);
	$sk	= 0;
	foreach $sn (sort numerically keys %mrna) {
		$sflag2	= 0;
		foreach $sn1 (sort numerically keys %mrna) {
			if ($sn1 <= $sn) {
				next;
			}
			$sflag	= -1;
			if ($mrna{$sn}{"L"}{"start"} <= $mrna{$sn1}{"L"}{"start"} && $mrna{$sn1}{"L"}{"start"} <= $mrna{$sn}{"L"}{"end"} && $mrna{$sn}{"L"}{"end"} <= $mrna{$sn1}{"L"}{"end"}) {
	#			print STDERR "site1: sn=",$mrna{$sn}{"L"}{"site"},",sn1=",$mrna{$sn1}{"L"}{"site"};
				$mrna{$sn}{"L"}{"end"}	= $mrna{$sn1}{"L"}{"end"};
#				$mrna{$sn}{"L"}{"hit"} = $mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"} ? $mrna{$sn1}{"L"}{"hit"} : $mrna{$sn}{"L"}{"hit"};
#				$mrna{$sn}{"L"}{"num"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"num"} : $mrna{$sn}{"L"}{"num"};
				$sflag	= $sn;	#	print STDERR "\tnew sn end=",$mrna{$sn}{"L"}{"end"},"\n";
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"L"}{"start"} >= $mrna{$sn1}{"L"}{"start"} && $mrna{$sn1}{"L"}{"end"} >= $mrna{$sn}{"L"}{"start"} && $mrna{$sn}{"L"}{"end"} >= $mrna{$sn1}{"L"}{"end"}) {
	#			print STDERR "site2: sn=",$mrna{$sn}{"L"}{"site"},",sn1=",$mrna{$sn1}{"L"}{"site"};
				$mrna{$sn}{"L"}{"start"}	= $mrna{$sn1}{"L"}{"start"};
#				$mrna{$sn}{"L"}{"hit"} = $mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"} ? $mrna{$sn1}{"L"}{"hit"} : $mrna{$sn}{"L"}{"hit"};
#				$mrna{$sn}{"L"}{"num"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"num"} : $mrna{$sn}{"L"}{"num"};
				$sflag	= $sn;	#	print STDERR "\tnew sn start=",$mrna{$sn}{"L"}{"start"},"\n";
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"L"}{"start"} <= $mrna{$sn1}{"L"}{"start"} && $mrna{$sn1}{"L"}{"end"} <= $mrna{$sn}{"L"}{"end"}) {
#				$mrna{$sn}{"L"}{"hit"} = $mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"} ? $mrna{$sn1}{"L"}{"hit"} : $mrna{$sn}{"L"}{"hit"};
#				$mrna{$sn}{"L"}{"num"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"num"} : $mrna{$sn}{"L"}{"num"};
	#			print STDERR "site3: sn=",$mrna{$sn}{"L"}{"site"},",sn1=",$mrna{$sn1}{"L"}{"site"},"\n";
				$sflag	= $sn;
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"L"}{"start"} >= $mrna{$sn1}{"L"}{"start"} && $mrna{$sn1}{"L"}{"end"} >= $mrna{$sn}{"L"}{"end"}) {
	#			print STDERR "site4: sn=",$mrna{$sn}{"L"}{"site"},",sn1=",$mrna{$sn1}{"L"}{"site"};
				$mrna{$sn}{"L"}{"end"}	= $mrna{$sn1}{"L"}{"end"};
				$mrna{$sn}{"L"}{"start"}	= $mrna{$sn1}{"L"}{"start"};
#				$mrna{$sn}{"L"}{"hit"} = $mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"} ? $mrna{$sn1}{"L"}{"hit"} : $mrna{$sn}{"L"}{"hit"};
#				$mrna{$sn}{"L"}{"num"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"num"} : $mrna{$sn}{"L"}{"num"};
				$sflag	= $sn;	#	print STDERR "\tnew sn start=",$mrna{$sn}{"L"}{"start"},"\tnew sn end=",$mrna{$sn}{"L"}{"end"},"\n";
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			}
			if ($sflag	!= -1) {
				if ($mrna{$sn}{"R"}{"start"} > $mrna{$sn1}{"R"}{"start"}) {
					$mrna{$sn}{"R"}{"start"}	= $mrna{$sn1}{"R"}{"start"};
				} elsif ($mrna{$sn}{"R"}{"end"} < $mrna{$sn1}{"R"}{"end"}) {
					$mrna{$sn}{"R"}{"end"}	= $mrna{$sn1}{"R"}{"end"};
				}
				if ($mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"}) {
					$mrna{$sn}{"L"}{"hit"} =  $mrna{$sn1}{"L"}{"hit"};
				}
				if ($mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"}) {
					$mrna{$sn}{"L"}{"num"} = $mrna{$sn1}{"L"}{"num"};
					$mrna{$sn}{"L"}{"site"} =  $mrna{$sn1}{"L"}{"site"};
				}
#				print STDERR "Lsite: sn=",$mrna{$sn}{"L"}{"site"},"\tsn1=",$mrna{$sn1}{"L"}{"site"},"\tLnum:sn=",$mrna{$sn}{"L"}{"num"},"\tsn1=",$mrna{$sn1}{"L"}{"num"},"\n";
				$mrna{$sn}{"R"}{"hit"} = $mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"} ? $mrna{$sn1}{"R"}{"hit"} : $mrna{$sn}{"R"}{"hit"};
				$mrna{$sn}{"R"}{"num"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"R"}{"num"} : $mrna{$sn}{"R"}{"num"};
#				print STDERR $mrna{$sn}{"L"}{"start"},"..",$mrna{$sn}{"L"}{"end"},"; L del: ",$mrna{$sn1}{"L"}{"start"},"..",$mrna{$sn1}{"L"}{"end"},"\t";
				$mrna{$sn1} = ();	delete($mrna{$sn1});	$sk++;
#				print STDERR "hit:L",$mrna{$sn}{"L"}{"hit"},",R",$mrna{$sn}{"R"}{"hit"},"\tnum:L",$mrna{$sn}{"L"}{"num"},",R",$mrna{$sn}{"R"}{"num"},"\n";
				$sflag2	= 1;	last;
			}
			$sflag	= -1;
			if ($mrna{$sn}{"R"}{"start"} <= $mrna{$sn1}{"R"}{"start"} && $mrna{$sn1}{"R"}{"start"} <= $mrna{$sn}{"R"}{"end"} && $mrna{$sn}{"R"}{"end"} <= $mrna{$sn1}{"R"}{"end"}) {
				$mrna{$sn}{"R"}{"end"}	= $mrna{$sn1}{"R"}{"end"};
#				$mrna{$sn}{"R"}{"hit"} = $mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"} ? $mrna{$sn1}{"R"}{"hit"} : $mrna{$sn}{"R"}{"hit"};
#				$mrna{$sn}{"R"}{"num"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"R"}{"num"} : $mrna{$sn}{"R"}{"num"};
				$sflag	= $sn;
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"R"}{"start"} >= $mrna{$sn1}{"R"}{"start"} && $mrna{$sn1}{"R"}{"end"} >= $mrna{$sn}{"R"}{"start"} && $mrna{$sn}{"R"}{"end"} >= $mrna{$sn1}{"R"}{"end"}) {
				$mrna{$sn}{"R"}{"start"}	= $mrna{$sn1}{"R"}{"start"};
#				$mrna{$sn}{"R"}{"hit"} = $mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"} ? $mrna{$sn1}{"R"}{"hit"} : $mrna{$sn}{"R"}{"hit"};
#				$mrna{$sn}{"R"}{"num"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"R"}{"num"} : $mrna{$sn}{"R"}{"num"};
				$sflag	= $sn;
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"R"}{"start"} <= $mrna{$sn1}{"R"}{"start"} && $mrna{$sn1}{"R"}{"end"} <= $mrna{$sn}{"R"}{"end"}) {
#				$mrna{$sn}{"R"}{"hit"} = $mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"} ? $mrna{$sn1}{"R"}{"hit"} : $mrna{$sn}{"R"}{"hit"};
#				$mrna{$sn}{"R"}{"num"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"R"}{"num"} : $mrna{$sn}{"R"}{"num"};
				$sflag	= $sn;
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			} elsif ($mrna{$sn}{"R"}{"start"} >= $mrna{$sn1}{"R"}{"start"} && $mrna{$sn1}{"R"}{"end"} >= $mrna{$sn}{"R"}{"end"}) {
				$mrna{$sn}{"R"}{"end"}	= $mrna{$sn1}{"R"}{"end"};
				$mrna{$sn}{"R"}{"start"}	= $mrna{$sn1}{"R"}{"start"};
#				$mrna{$sn}{"R"}{"hit"} = $mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"} ? $mrna{$sn1}{"R"}{"hit"} : $mrna{$sn}{"R"}{"hit"};
#				$mrna{$sn}{"R"}{"num"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"R"}{"num"} : $mrna{$sn}{"R"}{"num"};
				$sflag	= $sn;
#				$mrna{$sn}{"L"}{"site"} = $mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"} ? $mrna{$sn1}{"L"}{"site"} : $mrna{$sn}{"L"}{"site"};
			}
		#	print STDERR "Rsite: sn=",$mrna{$sn}{"L"}{"site"},"\tsn1=",$mrna{$sn}{"L"}{"site"},"\n";
			if ($sflag	!= -1) {
				if ($mrna{$sn}{"L"}{"start"} > $mrna{$sn1}{"L"}{"start"}) {
					$mrna{$sn}{"L"}{"start"}	= $mrna{$sn1}{"L"}{"start"};
				} elsif ($mrna{$sn}{"L"}{"end"} < $mrna{$sn1}{"L"}{"end"}) {
					$mrna{$sn}{"L"}{"end"}	= $mrna{$sn1}{"L"}{"end"};
				}
				if ($mrna{$sn}{"R"}{"hit"} < $mrna{$sn1}{"R"}{"hit"}) {
					$mrna{$sn}{"R"}{"hit"} =  $mrna{$sn1}{"R"}{"hit"};
				}
				if ($mrna{$sn}{"R"}{"num"} < $mrna{$sn1}{"R"}{"num"}) {
					$mrna{$sn}{"R"}{"num"} =  $mrna{$sn1}{"R"}{"num"};
					$mrna{$sn}{"L"}{"site"} = $mrna{$sn1}{"L"}{"site"};
				}
#				print STDERR "Rsite: sn=",$mrna{$sn}{"L"}{"site"},"\tsn1=",$mrna{$sn1}{"L"}{"site"},"\tRnum:sn=",$mrna{$sn}{"R"}{"num"},"\tsn1=",$mrna{$sn1}{"R"}{"num"},"\n";
				$mrna{$sn}{"L"}{"hit"} = $mrna{$sn}{"L"}{"hit"} < $mrna{$sn1}{"L"}{"hit"} ? $mrna{$sn1}{"L"}{"hit"} : $mrna{$sn}{"L"}{"hit"};
				$mrna{$sn}{"L"}{"num"} = $mrna{$sn}{"L"}{"num"} < $mrna{$sn1}{"L"}{"num"} ? $mrna{$sn1}{"L"}{"num"} : $mrna{$sn}{"L"}{"num"};
#				print STDERR $mrna{$sn}{"R"}{"start"},"..",$mrna{$sn}{"R"}{"end"},"; R del: ",$mrna{$sn1}{"R"}{"start"},"..",$mrna{$sn1}{"R"}{"end"},"\n";
				$mrna{$sn1} = ();	delete($mrna{$sn1});	$sk++;
#				print STDERR "hit:L",$mrna{$sn}{"L"}{"hit"},",R",$mrna{$sn}{"R"}{"hit"},"\tnum:L",$mrna{$sn}{"L"}{"num"},",R",$mrna{$sn}{"R"}{"num"},"\n";
				$sflag2	= 1;	last;
			}
		}
		if ($sflag2	== 1) {
			last;
		}
	}
	if ($sk	> 0) {
#		print STDERR "there are $sk link\n";
		sub_link_max();
	}
	return $sk;
}
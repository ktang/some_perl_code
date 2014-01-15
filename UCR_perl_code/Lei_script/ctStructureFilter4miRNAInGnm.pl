#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			ctStructureFilter
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.09.10
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.11.06
# Description:	in ct file of RNA second structure, to choose the sequence in which the stem length > $mismatch.
#**************************
# Version: 1.1		check the free energy, set the FE threshold: $maxFree
# Version: 2.0		bulges (size: one or two bases; number one or less)
# Version: 2.1		fix the output, additional output the stem seqeuences.
# Version: 2.2		deal with the multi-structure of one sequence (hybrid-ss-min --mfold=*,*,5 will output 5 max optimal structures for one sequence)
#**************************
#refer to bracStructureFilter.pl

# e-mail:highlei@gmail.com

my $version="2.2";
print STDERR ("\n============================================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:m:s:w:c:f:b:n:o:v:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $mismatch	= (defined $opt_m) ? $opt_m : 5;
my $slideLen	= (defined $opt_s) ? $opt_s : 1;
my $winLen		= (defined $opt_w) ? $opt_w : 21;
my $colInct		= (defined $opt_c) ? $opt_c : 4;	## the pair column
my $maxFree		= (defined $opt_f) ? $opt_f : -35;
my $bulges_size	= (defined $opt_b) ? $opt_b : 2;
my $bulges_num	= (defined $opt_n) ? $opt_n : 1;
my $outfile		= (defined $opt_o) ? $opt_o : "";
my $overlap		= (defined $opt_v) ? $opt_v : 0.5;

if ($opt_h){# || $aafile eq "" || $ntfile eq "") { #|| $batchfile eq "" || $output eq "") {
	usage();
}
sub numerically{$a<=>$b};
#sub sub_slideWindow;
#sub sub_numOfDiff;

use FileHandle;
use strict;



my ($i,$j,$k,$num,$len,$ntlen,$k1,$k2,$k3,$m,$n,$file,$line,$line1,$in,$match,$omatch,$a,$b,$end);
my (@buf,@tmp,@brac,@tmp1);
my (%aa,%nt,%ct,%mrna);
my $key		="";
my $key2	="";
my $flag	=0;
my ($opNum,$fore,$back,$gnmNum,$opRNA,$opRNA0);
my ($ovlp_s,$ovlp_e);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
#	print "retro-translate amino-acid sequences (alignment) to nucleotide sequence\n";
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %40s : %-25s","input RNAfold result file name (ct)",$infile);#%45s
	printf STDERR ("\n m  %40s : %-25s","number of mismatch in stem",$mismatch);
	printf STDERR ("\n s  %40s : %-25s","length of slide each time",$slideLen);
	printf STDERR ("\n w  %40s : %-25s","length of window",$winLen);
	printf STDERR ("\n f  %40s : %-25s","allowable max free energy",$maxFree);
	printf STDERR ("\n b  %40s : %-25s","allowable max asymmetric bulge size",$bulges_size);
	printf STDERR ("\n n  %40s : %-25s","allowable max asymmetric bulges num",$bulges_num);
	printf STDERR ("\n o  %40s : %-25s","output the possible mature miRNA sequences",$outfile);
	printf STDERR ("\n x  %40s !","exit the program");
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input RNAfold result file name (ct):\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "m") {print STDERR "please input the number of mismatch in stem\n"; $mismatch	= <STDIN>;$mismatch	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input length of slide each time\n"; $slideLen	= <STDIN>;$slideLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "w") {print STDERR "please input length of window\n"; $winLen	= <STDIN>;$winLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input allowable max free energy\n"; $maxFree	= <STDIN>;$maxFree	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "b") {print STDERR "please input allowable max asymmetric bulge size\n"; $bulges_size	= <STDIN>;$bulges_size	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "n") {print STDERR "please input allowable max asymmetric bulge number\n"; $bulges_num	= <STDIN>;$bulges_num	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "o") {print STDERR "please input file name to store the possible mature miRNA sequences\n"; $outfile	= <STDIN>;$outfile	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {sub_end_program();;exit(0);}
}


#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";
############################################ read ct files ######################################################

$i      = 0;	$n	= 0;
$file = new FileHandle ("$infile") || die;
while (<$file>) {
	$_=~s/^[\s|\t]+//g;
	$_=~s/[\s|\t|\r|\n]+$//g;
        if ($_ =~/\=/) {
			$j      =$_;    splice(@buf,0);
			@buf    = split(/\t+/,$j);
			$k	=	$buf[2];		## name of sequence
#				print STDERR "name=$k,    ";
			if (!exists($ct{$k})) {
				$ct{$k}{"num"}		= 1;	$n	= 1;
			} else {
				$ct{$k}{"num"}++;	$n	= $ct{$k}{"num"};	#print STDERR "n=$n,num=",$ct{$k}{"num"},"\n";
			}
			$ct{$k}{$n}{"title"}	= $j;
			$j	=~ /\=\s?(\S+)/;
			$ct{$k}{$n}{"FE"}	= $1;
			$i++;	$k2 = 0;	$i%10000 != 0 || print STDERR "$i, ";
        } else {
			$j	= $_;
			splice(@buf,0);
			@buf    = split(/\s+/,$j);
#			print STDERR "buf=@buf,    ";
			for ($k1 = 0; $k1 < @buf;$k1++) {
				$ct{$k}{$n}{$k2}{$k1}	= $buf[$k1];
			}
#			print STDERR $ct{$k}{$n}{$k2}{4},", ";
			$k2++;
		}
}
close $file || die("Wrong!");

print STDERR "\nThere are ",$i," sequences in file: $infile\n";
############################################ check the sequence ######################################################
$i	= 0;
foreach $key (sort keys %ct) {
#	print STDERR "num=",$ct{$key}{"num"},"\n";
	foreach $k (sort numerically keys %{$ct{$key}}) {
		if ($k eq "num") {
			next;
		}
		if ($ct{$key}{$k}{"FE"} > $maxFree) {
	#		print STDERR "del: $key,$k,",$ct{$key}{$k}{"title"},"\n";
			delete($ct{$key}{$k});
			$ct{$key}{"num"}--;
			$i++;
		}
	}
}
print STDERR "\nThere are ",$i," sequences which are deleted due to free energy > $maxFree\n";
$i	= 0;
foreach $key (sort keys %ct) {
	$i++;
}
print STDERR "There are ",$i," sequences remained after check free energy\n\n";

if ($outfile ne "") {
	open(OUTF,">$outfile") || die("Cannot open file: $outfile\n!");
}

$i	= 0;	$j++;
foreach $key (sort keys %ct) {
	$line	= "";	$opRNA0	= "";
	foreach $key2 (sort numerically keys %{$ct{$key}}) {
		if ($key2 eq "num") {
			next;
		}
		$ct{$key}{$key2}{"title"}	=~ /(\d+)\,(\d+)\s*$/;
		$ovlp_s	= $1;	$ovlp_e	= $2;	#print STDERR $ct{$key}{$key2}{"title"},"\ts=$ovlp_s\te=$ovlp_e\n";
		splice(@buf,0);				## the pair column
		foreach $k (sort numerically keys %{$ct{$key}{$key2}}) {
			if ($k eq "title" || $k eq "FE") {
				next;
			}
			$buf[$k]	= $ct{$key}{$key2}{$k}{$colInct};
	#		print STDERR "$buf[$k], ";
		}
	#	print STDERR "buf=@buf,    ";
		$line1	= "";	#$opRNA	= "";
		for ($k1 = 0;$k1 <= @buf - $winLen ;$k1+=$slideLen) {
			splice(@tmp,0);	splice(@tmp1,0);	$opRNA	= "";
			for ($k2 = $k1;$k2 < $k1+$winLen ;$k2++) {
				$tmp[$k2 - $k1]	= $buf[$k2];
			}
	#		@tmp = substr(@buf,$k1,$winLen);
	#		print STDERR "tmp=@tmp,    ";
			#----------------------------------------------------mismatch---------------------------------
			$m	= 0;
			for ($k2 = 0;$k2 < $winLen ;$k2++) {
				if ($tmp[$k2] == 0) {
					$m++;
				}
				$tmp1[$k2]	= $tmp[$k2];
			}
			if ($m > $mismatch) {						# m = mismatch
				next;
			}
			#----------------------------------------------------order---------------------------------
	#		print STDERR "mismatch=$m, ";
			for ($k2 = 0;$k2 < @tmp1 ;$k2++) {
				if ($tmp1[$k2] == 0) {
					splice(@tmp1,$k2,1);
					$k2--;
				}
			}
	#		$m	= @tmp;
	#		print STDERR "len=$m,tmp=@tmp;    ";

			$m	= 0;	$n	= 0;						# n,m = from large to small or reverse
			for ($k2 = 1;$k2 < @tmp1 ;$k2++) {
				if ($tmp1[$k2] > $tmp1[$k2-1]) {
					$m	= 1;							# m = from small to large
				} elsif ($tmp1[$k2] < $tmp1[$k2 -1]) {
					$n	= 1;							# n = from large to small
				}
			}
			$flag0	= 0;
			if ($m > 0 && $n > 0) {
				next;
			} elsif ($m > 0) {
				$flag0	= 1;							# m = from small to large
			} elsif ($n	> 0) {
				$flag0	= 2;							# n = from large to small
			}
			#----------------------------------------------------mismatch---------------------------------
			@tmp1	= sort numerically @tmp1;
			$m	= $tmp1[@tmp1-1];	$n	= $tmp1[0];		# m, n = max and min in array @tmp1
			if ( $m - $n < $winLen - $mismatch) {
				next;
			}
			if ($m - $n > $winLen + $mismatch) {
				next;
			}
			if ($m <= $k1 || $n >= $k1+$winLen) {
			} else {
				next;
			}
			#----------------------------------------------------bulges---------------------------------
			if ($flag0	== 1) {
				$m	= 1;	$flag	= 0;	$n	= 0;
				for ($k2 = 0;$k2 < $winLen ;$k2++) {
					if ($tmp[$k2] != 0) {
						last;
					}
				}
				for ($k2++;$k2 < $winLen ;$k2++) {
					if ($tmp[$k2] == 0) {
						$m++; next;
					}
		#			if ($tmp[$k2]-$tmp[$k2-$m] > $bulges_size+$m) {
		#				$flag	= 1;	last;
		#			} elsif ($tmp[$k2]-$tmp[$k2-$m] > $m) {
		#				$n++;
		#			} elsif ($tmp[$k2]-$tmp[$k2-$m] < $m-$bulges_size) {
		#				$flag	= 1;	last;
		#			} elsif ($tmp[$k2]-$tmp[$k2-$m] < $m) {
		#				$n++;
		#			}	
					if ($tmp[$k2]-$tmp[$k2-$m] > $m) {
						if ($tmp[$k2]-$tmp[$k2-$m] > $bulges_size) {
							$flag	= 1;	last;
						}
						$n++;
					} elsif ($tmp[$k2]-$tmp[$k2-$m] < $m) {
						if ($tmp[$k2]-$tmp[$k2-$m] < $m-$bulges_size) {
							$flag	= 1;	last;
						} 
						$n++;
					} 
					$m	= 1;
				}
				if ($flag > 0 || $n > $bulges_num) {
					next;
				}

			} elsif ($flag0	== 2) {							# 2 = from large to small
				$m	= 1;	$flag	= 0;	$n	= 0;
				for ($k2 = 0;$k2 < $winLen ;$k2++) {
					if ($tmp[$k2] != 0) {
						last;
					}
				}
				for ($k2++;$k2 < $winLen ;$k2++) {
					if ($tmp[$k2] == 0) {
						$m++; next;
					}	
					if ($tmp[$k2-$m]-$tmp[$k2] > $m) {
						if ($tmp[$k2-$m]-$tmp[$k2] > $bulges_size) {
							$flag	= 1;	last;
						}
						$n++;
					} elsif ($tmp[$k2-$m]-$tmp[$k2] < $m) {
						if ($tmp[$k2-$m]-$tmp[$k2] < $m-$bulges_size) {
							$flag	= 1;	last;
						} 
						$n++;
					} 
					$m	= 1;
				}
				if ($flag > 0 || $n > $bulges_num) {
					next;
				}
			} else {
				print STDERR "m=$m,flag=$flag is wrong!\n";
			}
			#----------------------------------------------------bulges---------------------------------
			$m	= sub_overlap($k1+1,$k1+$winLen,$tmp1[0],$tmp1[@tmp1-1],$ovlp_s,$ovlp_e,$overlap);
			if ($m == 0) {
				next;
			}
			$line1	.= ($k1+1). ".." . ($k1+$winLen);
			$opRNA	.= 	">".$key."_" . ($k1+1). ".." . ($k1+$winLen) . "\n";
			for ($m = $k1;$m < $k1+$winLen ;$m++) {
				$opRNA	.= $ct{$key}{$key2}{$m}{1};
			}
			$opRNA	.= "\n";
			if (!exists($mrna{$opRNA})) {
				$mrna{$opRNA}	= 1;
				$opRNA0	.= $opRNA;
			}
			if ($flag0	== 1) {
				$line1	.= ",".$tmp1[0]."..".$tmp1[@tmp1-1]."; ";
			} elsif ($flag0	== 2) {
				$line1	.= ",".$tmp1[@tmp1-1]."..".$tmp1[0]."; ";
			} 
			$j++;	#last;
		}
		if ($line1 =~ /\,/) {
			$line	.= "\n".$key."\t" . $line1;
		}
	}
	if ($line	=~ /\,/) {
		print substr($line,1),"\n";	$i++;
		if ($outfile ne "" && $opRNA0 ne "") {
			print OUTF	$opRNA0;
		}
	}
}
if ($outfile ne "") {
	close(OUTF);
}
print STDERR "There are $i,$j sequences,pieces in output\n";

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
#	print "retro-translate amino-acid sequences (alignment) to nucleotide sequences\n";
	print "Usage:\n	$0 \n";
	print "-i	input RNAfold result file name (ct)\n";
	print "			eg: RNAfold.brac\n";
	print "-m	input number of mismatch in stem\n";
	print "			eg: -m 4; default: 5\n";
	print "-s	input length of slide each time\n";
	print "			eg: -s 2; default: 1\n";
	print "-w	please input length of window\n";
	print "			eg: -w 22; default: 21\n";
	print "-f	please input allowable max free energy\n";
	print "			eg: -f -40; default: -35\n";
	print "-b	please input allowable max asymmetric bulge size\n";
	print "			eg: -b 1; default: 2\n";
	print "-n	please input allowable max asymmetric bulges number\n";
	print "			eg: -n 0; default: 1\n";
	print "-o	output the possible mature miRNA sequences\n";
	print "			eg: -o possible_miRNA.fa; default: donot output\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i RNAfold.ct -m 4\n";
	print "$0 -i RNAfold.ct -m 5 -s 1 -w 21\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 1\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 1 -o possible.miRNA.fa\n";
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
######################                  sub_overlap
############################################################################################################

sub sub_overlap
{
	my ($s,$e,$ps,$pe,$ovs,$ove,$opRatio)	= @_;
	my $vv=0;
	if ($s > $ovs && $s < $ove && $e > $ove) {
		$vv	= 1.0 * ($ove-$s+1)/($e-$s+1);
	} elsif ($s < $ovs && $ovs < $e && $e < $ove) {
		$vv	= 1.0 * ($e-$ovs+1)/($e-$s+1);
	} elsif ($s < $ovs && $ove < $e) {
		$vv	= 1;
	} elsif ($s > $ovs && $ove > $e) {
		$vv	= 1;
	}
	if ($vv >= $opRatio) {
		return 1;
	}
	$s	= $ps;	$e	= $pe;
	if ($s > $ovs && $s < $ove && $e > $ove) {
		$vv	= 1.0 * ($ove-$s+1)/($e-$s+1);
	} elsif ($s < $ovs && $ovs < $e && $e < $ove) {
		$vv	= 1.0 * ($e-$ovs+1)/($e-$s+1);
	} elsif ($s < $ovs && $ove < $e) {
		$vv	= 1;
	} elsif ($s > $ovs && $ove > $e) {
		$vv	= 1;
	}
	if ($vv >= $opRatio) {
		return 1;
	}
	return 0;
}

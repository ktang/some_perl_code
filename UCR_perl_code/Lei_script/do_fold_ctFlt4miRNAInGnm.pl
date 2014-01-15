#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			do_fold_ctFlt4miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2009.11.18
# Description:	read fasta file, then do unafold, then do ctStrucFilter
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
getopts("hi:u:m:s:w:c:f:b:n:r:v:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 0;
my $infile		= $opt_i;
my $fold_set	= (defined $opt_u) ? $opt_u : "" ; # --mfold=*,*,5

my $mismatch	= (defined $opt_m) ? $opt_m : 5;
my $slideLen	= (defined $opt_s) ? $opt_s : 1;
my $winLen		= (defined $opt_w) ? $opt_w : 21;
my $colInct		= (defined $opt_c) ? $opt_c : 4;	## the pair column
my $maxFree		= (defined $opt_f) ? $opt_f : -35;
my $bulges_size	= (defined $opt_b) ? $opt_b : 2;
my $bulges_num	= (defined $opt_n) ? $opt_n : 1;
my $outRNA		= (defined $opt_r) ? $opt_r : "";
#my $overlap		= (defined $opt_v) ? $opt_v : 0.5;

if ($opt_h){
	usage();
}

my $unafold		= "hybrid-ss-min -s DAT";

sub numerically{$a<=>$b};

use FileHandle;
use strict;



my ($i,$j,$k,$num,$len,$ntlen,$k1,$k2,$k3,$m,$n,$file,$line,$line1,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@genome,@gnmName,@gnmLen);
my (%ct,%mrna);
my $key="";
my ($ii,$jj);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %40s : %-25s","input fasta file name",$infile);#%45s
	printf STDERR ("\n u  %45s : %-25s","the parameter set for unafold",$fold_set);
	printf STDERR ("\n m  %40s : %-25s","number of mismatch in stem",$mismatch);
	printf STDERR ("\n s  %40s : %-25s","length of slide each time",$slideLen);
	printf STDERR ("\n w  %40s : %-25s","length of window",$winLen);
	printf STDERR ("\n f  %40s : %-25s","allowable max free energy",$maxFree);
	printf STDERR ("\n b  %40s : %-25s","allowable max asymmetric bulge size",$bulges_size);
	printf STDERR ("\n n  %40s : %-25s","allowable max asymmetric bulges num",$bulges_num);
	printf STDERR ("\n r  %40s : %-25s","output the possible mature miRNA sequences",$outRNA);
	printf STDERR ("\n x  %40s !","exit the program");

	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input fasta file name:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "u") {print STDERR "please input the parameter set for unafold\n";$fold_set	= <STDIN>;$fold_set	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "m") {print STDERR "please input the number of mismatch in stem\n"; $mismatch	= <STDIN>;$mismatch	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input length of slide each time\n"; $slideLen	= <STDIN>;$slideLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "w") {print STDERR "please input length of window\n"; $winLen	= <STDIN>;$winLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input allowable max free energy\n"; $maxFree	= <STDIN>;$maxFree	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "b") {print STDERR "please input allowable max asymmetric bulge size\n"; $bulges_size	= <STDIN>;$bulges_size	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "n") {print STDERR "please input allowable max asymmetric bulge number\n"; $bulges_num	= <STDIN>;$bulges_num	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "r") {print STDERR "please input file name to store the possible mature miRNA sequences\n"; $outRNA	= <STDIN>;$outRNA	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {sub_end_program();;exit(0);}
}

#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";

############################################ read file ######################################################

############################################ read genome files ######################################################
$i	= -1;
$file = new FileHandle ("$infile") || die;
while (<$file>) {
	if ($_ =~/^>(\S+)/) {
			$i++;
			$genome[$i]	= "";
			$gnmName[$i]	= $_;	#$k1++;
	} else {
		$genome[$i] .= $_;
	}
}
close $file || die("Wrong!");
print STDERR "There are ",$i+1," sequences in file $infile\n";
############################################ cut genome and do fold ######################################################

open(CTFLT, ">$infile.ctFlt.txt") || die("Can not open file: $infile.ctFlt.txt\n");
open(STEMF, ">$infile.ctFlt.stem") || die("Can not open file: $infile.ctFlt.stem\n");
open(CTF, ">$infile.ctFlt.ct") || die("Can not open file: $infile.ctFlt.ct\n");

$ii	= 0;	$jj	= 0;
$j	= "__do_fold_$infile\_" . $start . ".fa";
for ($k2=0; $k2 < @genome;$k2++) {
	open(TMPF, ">$j") || die("Can not open file: $j\n");
	print TMPF $gnmName[$k2],$genome[$k2],"\n";
	close(TMPF);
	system("$unafold $fold_set $j > $j.unafold.out");
	%ct	= ();
	sub_readCTfile($j.".ct");
	$i	= 0;
	foreach $key (sort numerically keys %ct) {
		if ($key eq "num") {
			next;
		}
		if ($ct{$key}{"FE"} > $maxFree) {
			delete($ct{$key});
			$ct{"num"}--;
			$ii++;
		}
	}
#	print STDERR "\nThere are ",$i," sequences which are deleted due to free energy > $maxFree\n";
#	print STDERR "There are ",$ct{"num"}," sequences remained after check free energy\n\n";
	if ($ct{"num"}	== 0) {
		next;
	}
	($i,$a,$b,$file)	= sub_strucCheck();			# i=num, a=ctflt.txt, b=stem, file=ctfile.
	if ($i	== 0) {
		next;
	}
	$jj	+= $i;
	print	CTFLT	$a;		#,"\n";
	print	STEMF	$b;		#,"\n";
	print	CTF		$file;	#,"\n";
#	$i	= <STDIN>;
}
close(CTFLT) || die("Can not close file: $infile.ctFlt.txt\n");
close(STEMF) || die("Can not close file: $infile.stem\n");
close(CTF) || die("Can not close file: $infile.ct\n");
system("rm $j*");

print STDERR "\nThere are ",$ii," sequences which are deleted due to free energy > $maxFree.\n";
print STDERR "There are ",$jj," sequences remained after structure check.\n\n";




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
	print "-i	input RNAfold result file name (ct)\n";
	print "			eg: RNAfold.brac\n";
	print "-u	the parameter set for unafold\n";
	print "			eg: -u --mfold=*,*,5; default: $fold_set\n";
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
	print "-r	output the possible mature miRNA sequences\n";
	print "			eg: -o possible_miRNA.fa; default: donot output\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i RNAfold.ct -m 4\n";
	print "$0 -i RNAfold.ct -m 5 -s 1 -w 21\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 1\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 1 -r possible.miRNA.fa\n";
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
	printf STDERR ("do_fold_ctFlt4miRNAInGnm.pl execute time : %.2f s\n",$end-$start);
	print STDERR ("............................................................\n");
	exit(0);
	return;
}

############################################################################################################
######################                  sub_readCTfile
############################################################################################################

sub sub_readCTfile
{
	my ($rcname) = @_;
	my ($rcj,$rci,$rcm,$rcn,$rcfile);
	my (@rcbuf);

	$rcfile = new FileHandle ("$rcname") || die("Cannot open ct file: $rcname!\n");
	while (<$rcfile>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		splice(@rcbuf,0);
		if ($_ =~/\=/) {
			$rcj			= $_;
			@rcbuf			= split(/\t+/,$rcj);
			if (!exists($ct{"num"})) {
				$ct{"num"}		= 1;	$rcn	= 1;
			} else {
				$ct{"num"}++;	$rcn	= $ct{"num"};	#print STDERR "n=$n,num=",$ct{$k}{"num"},"\n";
			}
			$ct{$rcn}{"name"}		= $rcbuf[2];		## name of sequence
			$ct{$rcn}{"len"}		= $rcbuf[0];
			$ct{$rcn}{"title"}		= $rcj;		#	print STDERR $ct{$rcn}{"title"},"=title in read ct file\n";
			$rcj			=~ /\=\s?(\S+)/;
			$ct{$rcn}{"FE"}		= $1;
			$rci			= 0;
		} else {
			$rcj	= $_;
			@rcbuf    = split(/\s+/,$rcj);
			for ($rcm = 0; $rcm < @rcbuf;$rcm++) {
				$ct{$rcn}{$rci}{$rcm}	= $rcbuf[$rcm];
			}
			$rci++;
		}
	}
	close $rcfile || die("Wrong!");

	return $rcn;
}

############################################################################################################
######################                  sub_strucCheck
############################################################################################################

sub sub_strucCheck
{
	my ($sk,$sj,$skey,$sk1,$sk2,$sm,$sn,$sline,$sline1,$sflag,$sflag0,$opRNA,$opRNA0,$sct);	#,$ovlp_s,$ovlp_e
	my (@sbuf,@stmp,@stmp1);
	my %smrna=();
	my ($ovlp_s,$ovlp_e,$swinLen);

	$sline	= "";	$opRNA0	= "";	$sj	= 0;	$sct	= "";
	foreach $skey (sort numerically keys %ct) {
		if ($skey eq "num") {
			next;
		}
		$ct{$skey}{"title"}	=~ /(\d+)\,(\d+)\s*$/;
		$ovlp_s	= $1;	$ovlp_e	= $2;	#print STDERR $ct{$skey}{"title"},"\ts=$ovlp_s\te=$ovlp_e\n";
		$swinLen	= $ovlp_e - $ovlp_s + 1 > $winLen ? $winLen : $ovlp_e - $ovlp_s + 1;	#print STDERR "swinLen=$swinLen\n";
		splice(@sbuf,0);				## the pair column
		foreach $sk (sort numerically keys %{$ct{$skey}}) {
			if ($sk eq "title" || $sk eq "FE" || $sk eq "name" || $sk eq "len") {
				next;
			}
			$sbuf[$sk]	= $ct{$skey}{$sk}{$colInct};
		}
		$sline1	= "";	#$opRNA	= "";
		for ($sk1 = 0;$sk1 <= @sbuf - $swinLen ;$sk1+=$slideLen) {
			splice(@stmp,0);	splice(@stmp1,0);	$opRNA	= "";
			$sm	= 0;	$sn	= 0;
			for ($sk2 = $sk1;$sk2 < $sk1+$swinLen ;$sk2++) {
				$stmp[$sk2 - $sk1]	= $sbuf[$sk2];
				if ($sbuf[$sk2]!=0 && ($sbuf[$sk2] < $ovlp_s || $sbuf[$sk2] > $ovlp_e)) {
					$sm	= 1;
				}
				if ($sk2 < $ovlp_s-1 || $sk2 > $ovlp_e-1) {
					$sn	= 1;
				}
			}
			if ($sm == 1 && $sn	== 1) {
				next;
			}
			#----------------------------------------------------mismatch---------------------------------
			$sm	= 0;
			for ($sk2 = 0;$sk2 < $swinLen ;$sk2++) {
				if ($stmp[$sk2] == 0) {
					$sm++;
				}
				$stmp1[$sk2]	= $stmp[$sk2];
			}
			if ($sm > $mismatch) {						# m = mismatch
				next;
			}
			#----------------------------------------------------order---------------------------------
	#		print STDERR "mismatch=$sm, ";
			for ($sk2 = 0;$sk2 < @stmp1 ;$sk2++) {
				if ($stmp1[$sk2] == 0) {
					splice(@stmp1,$sk2,1);
					$sk2--;
				}
			}
			$sm	= 0;	$sn	= 0;						# n,m = from large to small or reverse
			for ($sk2 = 1;$sk2 < @stmp1 ;$sk2++) {
				if ($stmp1[$sk2] > $stmp1[$sk2-1]) {
					$sm	= 1;							# m = from small to large
				} elsif ($stmp1[$sk2] < $stmp1[$sk2 -1]) {
					$sn	= 1;							# n = from large to small
				}
			}
			$sflag0	= 0;
			if ($sm > 0 && $sn > 0) {
				next;
			} elsif ($sm > 0) {
				$sflag0	= 1;							# m = from small to large
			} elsif ($sn	> 0) {
				$sflag0	= 2;							# n = from large to small
			}
			#----------------------------------------------------mismatch---------------------------------
			@stmp1	= sort numerically @stmp1;
			$sm	= $stmp1[@stmp1-1];	$sn	= $stmp1[0];		# m, n = max and min in array @stmp1
			if ( $sm - $sn < $swinLen - $mismatch) {
				next;
			}
			if ($sm - $sn > $swinLen + $mismatch) {
				next;
			}
			if ($sm <= $sk1 || $sn >= $sk1+$swinLen) {
			} else {
				next;
			}
			#----------------------------------------------------bulges---------------------------------
			if ($sflag0	== 1) {
				$sm	= 1;	$sflag	= 0;	$sn	= 0;
				for ($sk2 = 0;$sk2 < $swinLen ;$sk2++) {
					if ($stmp[$sk2] != 0) {
						last;
					}
				}
				for ($sk2++;$sk2 < $swinLen ;$sk2++) {
					if ($stmp[$sk2] == 0) {
						$sm++; next;
					}
					if ($stmp[$sk2]-$stmp[$sk2-$sm] > $sm) {
						if ($stmp[$sk2]-$stmp[$sk2-$sm] > $bulges_size) {
							$sflag	= 1;	last;
						}
						$sn++;
					} elsif ($stmp[$sk2]-$stmp[$sk2-$sm] < $sm) {
						if ($stmp[$sk2]-$stmp[$sk2-$sm] < $sm-$bulges_size) {
							$sflag	= 1;	last;
						} 
						$sn++;
					} 
					$sm	= 1;
				}
				if ($sflag > 0 || $sn > $bulges_num) {
					next;
				}

			} elsif ($sflag0	== 2) {							# 2 = from large to small
				$sm	= 1;	$sflag	= 0;	$sn	= 0;
				for ($sk2 = 0;$sk2 < $swinLen ;$sk2++) {
					if ($stmp[$sk2] != 0) {
						last;
					}
				}
				for ($sk2++;$sk2 < $swinLen ;$sk2++) {
					if ($stmp[$sk2] == 0) {
						$sm++; next;
					}	
					if ($stmp[$sk2-$sm]-$stmp[$sk2] > $sm) {
						if ($stmp[$sk2-$sm]-$stmp[$sk2] > $bulges_size) {
							$sflag	= 1;	last;
						}
						$sn++;
					} elsif ($stmp[$sk2-$sm]-$stmp[$sk2] < $sm) {
						if ($stmp[$sk2-$sm]-$stmp[$sk2] < $sm-$bulges_size) {
							$sflag	= 1;	last;
						} 
						$sn++;
					} 
					$sm	= 1;
				}
				if ($sflag > 0 || $sn > $bulges_num) {
					next;
				}
			} else {
				print STDERR "m=$sm,flag=$sflag is wrong!\n";
			}
			#----------------------------------------------------bulges---------------------------------
	#		$sm	= sub_overlap($sk1+1,$sk1+$swinLen,$stmp1[0],$stmp1[@stmp1-1],$ovlp_s,$ovlp_e,$overlap);
	#		if ($sm == 0) {
	#			next;
	#		}
			$sline1	.= ($sk1+1). ".." . ($sk1+$swinLen);
			$opRNA	.= 	">".$ct{$skey}{"name"}."_" . ($sk1+1). ".." . ($sk1+$swinLen) . "\n";
			for ($sm = $sk1;$sm < $sk1+$swinLen ;$sm++) {
				$opRNA	.= $ct{$skey}{$sm}{1};
			}
			$opRNA	.= "\n";
			if (!exists($smrna{$opRNA})) {
				$smrna{$opRNA}	= 1;
				$opRNA0	.= $opRNA;
			}
			if ($sflag0	== 1) {
				$sline1	.= ",".$stmp1[0]."..".$stmp1[@stmp1-1]."; ";
			} elsif ($sflag0	== 2) {
				$sline1	.= ",".$stmp1[@stmp1-1]."..".$stmp1[0]."; ";
			} 
			#last;
		}
		if ($sline1 =~ /\,/) {
			$sj++;
			$sline	.= "\n".$ct{$skey}{"name"}."\t" . $sline1;
			$sct	.= $ct{$skey}{"title"}."\n";	#print STDERR $ct{$skey}{"title"},"=title,ct=$sct----\n";
			foreach $sk (sort numerically keys %{$ct{$skey}}) {
				if ($sk eq "title" || $sk eq "FE" || $sk eq "name" || $sk eq "len") {
					next;
				}
				foreach $sm (sort numerically keys %{$ct{$skey}{$sk}}) {
					$sct	.= $ct{$skey}{$sk}{$sm}."\t";
				}
		#		$sct	.= "\n";
				$sct	= substr($sct,0,length($sct)-1)."\n";
			}
		#	print STDERR "ct=$sct\n";
		}
	}
#	print STDERR "ctflt=$sline\n";
	if ($sline	=~ /\,/) {
		$sn	= substr($sline,1)."\n";
	} else {
		$sn	= "";
	}
	return $sj,$sn,$opRNA0,$sct;
}

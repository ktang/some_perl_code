#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			pickupFromh2p4miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.03
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	20010.08.31
# Description:	read histogram2predict4miRNAInGnm.pl result and pick out the most possible miRNAs.
#**************************
# Version: 6.5s	plus/minus	consist with miRNAInGnm
# Version: 6.6s	add $fasta, use the copy number in small RNA database as the criteria.
#**************************
#refer to histogram2predict4miRNAInGnm.pl

# e-mail:highlei@gmail.com

my $version="6.6s";
print STDERR ("\n==========================| $0 start |==================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:f:c:n:r:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 1;
my $infile		= $opt_i;
my $h2pFasta	= $opt_f;
#my $minHitsNum	= (defined $opt_n) ? $opt_n : 5;
#my $minRatio	= (defined $opt_r) ? $opt_r : 0.75;
#my $winLen		= (defined $opt_w) ? $opt_w : 21;


if ($opt_h) { #|| $batchfile eq "" || $output eq "") {
	usage();
}
sub numerically{$a<=>$b};
#sub sub_slideWindow;
#sub sub_numOfDiff;
my (%ctFlt,%mrna,%pgnm);

use FileHandle;
use strict;



my ($i,$j,$k,$num,$len,$ntlen,$k1,$k2,$k3,$k4,$m,$n,$file,$line,$in,$match,$omatch,$a,$b,$end);
my (@buf,@bufo,@genome,@gnmName);
my (%op,%fasta);
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
	printf STDERR ("\n i  %40s : %-25s","input histogram2predict.pl result file name",$infile);#%45s
	printf STDERR ("\n f  %40s : %-25s","input h2p fasta file",$h2pFasta);#%45s
#	printf STDERR ("\n n  %40s : %-25s","input minnimum hits number",$minHitsNum);
#	printf STDERR ("\n r  %40s : %-25s","input minnimum ratio of duplex",$minRatio);
	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input histogram2predict.pl result file name:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input h2p fasta file:\n"; $h2pFasta	= <STDIN>;	$h2pFasta	=~s/[\s|\t|\r|\n]+$//g;}
#	elsif($yesorno eq "n") {print STDERR "please input minnimum hits number\n"; $minHitsNum	= <STDIN>;$minHitsNum	=~s/[\s|\t|\r|\n]+$//g;}
#	elsif($yesorno eq "r") {print STDERR "please input minimum ratio of duplex in segments:\n";$minRatio	= <STDIN>;$minRatio	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "x") {sub_end_program();;exit(0);}
}


#print STDERR "\npro=$pro\tparam=$param\tbachfile=$batchfile\top=$output\tdir=$dir\n";

############################################ read h2p files ######################################################
sub_readFastaFile($h2pFasta);

sub_pickup_from_h2p($infile);

foreach $k1 (sort keys %pgnm) {
	foreach $k2 (sort numerically keys %{$pgnm{$k1}} ) {
		foreach $k3 (sort keys %{$pgnm{$k1}{$k2}} ) {
			if (!exists($op{$pgnm{$k1}{$k2}{$k3}{"Line"}})) {
				$op{$pgnm{$k1}{$k2}{$k3}{"Line"}}	= 0;
	#			print STDERR "Line=",$pgnm{$k1}{$k2}{$k3}{"Line"};
			} 
		}
	}
}

foreach $i (sort keys %op) {
	print $i;
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
	print "-i	input histogram2predict.pl result file name\n";
	print "			eg: hist2pred.h2p\n";
	print "-f	input h2p fasta file\n";
	print "			eg: rice.h2p.fa\n";
#	print "-n	input minnimum hits number\n";
#	print "			eg: -n 30; default: $minHitsNum\n";
#	print "-r	input minimum ratio of duplex in segments\n";
#	print "			eg: -r 0.5; default: $minRatio\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i hist2pred.h2p -f rice.h2p.fa\n";# -m \"-k 3 -s DAT\"\n";
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
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}
############################################################################################################
######################                  sub_readFastaFile
############################################################################################################


sub sub_readFastaFile
{
	my ($rcname) = @_;
	my ($rcj,$rci,$rck,$rcm,$rcn,$rcfile,$rckk);
	my (@rcbuf);
	$rckk	= 0;
	$rcfile = new FileHandle ("$rcname") || die("Cannot open h2p fasta file: $rcname!\n");
	while (<$rcfile>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		if ($_ =~/^>(\S+).+\_(\d+)x\;\d+\,\d+\;-?\d+\,-?\d+$/) {
			$rckk++;#	$rci	= $1;
			$fasta{$1}	= $2;
		} else {
#			$fasta{$rci}	.= $_;
		}
	}
	close $rcfile || die("Wrong when close $rcname!");
	print STDERR "Load $rcname OK!\t $rckk\n";
	return $rckk;
}

############################################################################################################
######################                  sub_pickup_from_h2p
############################################################################################################

sub sub_pickup_from_h2p
{
	my ($pinfile)=@_;
	my ($pki,$pkj,$pk,$pk1,$pk2,$pk3,$pm,$pn,$pa,$pb,$pfile);
	my (@pbuf,@pst);
	my (%phash);

	$pk1	= 0;	$pkj	= -1;
	open (PINF,"$pinfile") || die("Cannot open file: $pinfile\n");
	$pkj	= 0;	#	print STDERR "$pinfile\n";
	while(<PINF>)
	{
		$pkj++;
		if ($_ =~/^([^\_]+)\_(\d+)\_(\d+)\t(\d+)\;(\S+)\t(\S+)\t(\S+)\s(\S+)/) {
			my $pchr = $1;	my $pstart = $2;	my $plen = $3;	my $pnumofseg = $4;
			my $pseg = $5;	#	my $pleft = $6;	my $pright = $7;
			my $pcov = $8;	my	$pline = $_;
			$pk2 = $2+$3;	$pki = "$1\_$2\_$3";	splice(@pbuf,0);
			
			@pbuf    = split(/\;\s?/,$pseg);	#print STDERR "@pbuf\n";
			$pa	= 0;	$pb	= $pstart+$plen;
			for ($pk3 = 0;$pk3 < @pbuf ;$pk3++) {
				if ($pbuf[$pk3] =~/(\d+)\.\.(\d+)\,(\d+)\.\.(\d+)/) {
					$pst[0] = $1;	$pst[1] = $2;	$pst[2] = $3;	$pst[3] = $4;
					@pst = sort numerically @pst;
					$pa	= $pst[3] > $pa ? $pst[3] : $pa;
					$pb = $pst[0] < $pb ? $pst[0] : $pb;
				} elsif ($pbuf[$pk3] =~/(-\d+)\.\.(-\d+)\,(-\d+)\.\.(-\d+)/) {
				#	$pst[0] = $1+$plen;	$pst[1] = $2+$plen;	$pst[2] = $3+$plen;	$pst[3] = $4+$plen;
					$pst[0] = 1-$1;	$pst[1] = 1-$2;	$pst[2] = 1-$3;	$pst[3] = 1-$4;
					@pst = sort numerically @pst;
					$pa	= $pst[3] > $pa ? $pst[3] : $pa;
					$pb = $pst[0] < $pb ? $pst[0] : $pb;
				} else {
					print STDERR "match wrong in 101:$pbuf[$pk3]\n";
				}
			}
			$pa += $pstart;	$pb += $pstart;		#print STDERR "$pa,$pb\n";
			my $pf0	= 1;	my $pf1 = 1;	my $pf2 = 1;	my $pf3	= 0;
			for ($pk3 = $pb;$pk3 <= $pa ;$pk3++) {
				if (exists($pgnm{$pchr}{$pk3})) {
					foreach $pk (sort keys %{$pgnm{$pchr}{$pk3}}) {
						if ($fasta{$pk} > $fasta{$pki}) {
							$pf0	= -1; last;
						} elsif ($fasta{$pk} < $fasta{$pki}) {
							$pf0	= 0;
						} else {
							$pf0	*= 1;
							if ($pgnm{$pchr}{$pk3}{$pk}{"Cov"} > $pcov) {
								$pf1	= -1; last;
							} elsif ($pgnm{$pchr}{$pk3}{$pk}{"Cov"} < $pcov) {
								$pf1	= 0;
							} elsif ($pgnm{$pchr}{$pk3}{$pk}{"Cov"} == $pcov) {
								$pf1	*= 1;
								if ($pgnm{$pchr}{$pk3}{$pk}{"Len"} > $plen) {
									$pf2	= 0; #print STDERR $pgnm{$pchr}{$pk3}{$pk}{"Len"},"=len,plen=$plen\n";
								} elsif ($pgnm{$pchr}{$pk3}{$pk}{"Len"} < $plen) {
									$pf2	= -1;	last;
								}elsif ($pgnm{$pchr}{$pk3}{$pk}{"Len"} == $plen) {
									$pf2	*= 1;
								}
							}
						}
					}
					if ($pf0	== -1 || $pf1 == -1 || $pf2 == -1) {
						last;
					}
					$pf3	= 1;
				} 
			}
			if ($pf0	== -1 || $pf1 == -1 || $pf2 == -1) {#	print STDERR "1\tpf1=$pf1,pf2=$pf2\n";
				next;
			}
			if ($pf3	== 0) {	# || $pf1 == 1 && $pf2 == 1
				#	print STDERR "2\tpf1=$pf1,pf2=$pf2,pf3=$pf3\n";
				for ($pk3 = $pb;$pk3 <= $pa ;$pk3++) {
					$pgnm{$pchr}{$pk3}{$pki}{"Len"}		= $plen;
					$pgnm{$pchr}{$pk3}{$pki}{"Segnum"}	= $pnumofseg;
					$pgnm{$pchr}{$pk3}{$pki}{"Line"}	= $pline;
					$pgnm{$pchr}{$pk3}{$pki}{"Cov"}		= $pcov;
					$pgnm{$pchr}{$pk3}{$pki}{"Start"}	= $pb;
					$pgnm{$pchr}{$pk3}{$pki}{"End"}		= $pa;
				}
				$pk1++;
				next;
			} 
			if ($pf0*$pf1*$pf2 == 1) {
				next;
			}
#			if ($pf1	== 0 || $pf1 == 1 && $pf2 != 1) {#	print STDERR "3\tpf1=$pf1,pf2=$pf2\n";
				$pm	= 0;	$pn	= $pa;
				for ($pk3 = $pb;$pk3 <= $pa ;$pk3++) {
					foreach $pk (sort keys %{$pgnm{$pchr}{$pk3}} ) {
						$pm = $pm > $pgnm{$pchr}{$pk3}{$pk}{"End"}   ? $pm : $pgnm{$pchr}{$pk3}{$pk}{"End"};
						$pn = $pn < $pgnm{$pchr}{$pk3}{$pk}{"Start"} ? $pn : $pgnm{$pchr}{$pk3}{$pk}{"Start"};
					}
				}
				$pf1 = 0;	$pf2 = 0;
				for ($pk3 = $pn; $pk3 <= $pm ; $pk3++) {
					if (exists($pgnm{$pchr}{$pk3})) {
						delete($pgnm{$pchr}{$pk3});
						if ($pf1 == 0) {
							$pf2++;	$pf1 = 1;
						}
					} else {
						if ($pf1 == 1) {
							$pf1 = 0;
						}
					}
				}
				$pk1	= $pk1 + 1 - $pf2;
				for ($pk3 = $pb;$pk3 <= $pa ;$pk3++) {
					$pgnm{$pchr}{$pk3}{$pki}{"Len"}		= $plen;
					$pgnm{$pchr}{$pk3}{$pki}{"Segnum"}	= $pnumofseg;
					$pgnm{$pchr}{$pk3}{$pki}{"Line"}	= $pline;
					$pgnm{$pchr}{$pk3}{$pki}{"Cov"}		= $pcov;
					$pgnm{$pchr}{$pk3}{$pki}{"Start"}	= $pb;
					$pgnm{$pchr}{$pk3}{$pki}{"End"}		= $pa;
				}
#				next;
#			}
			
		} else {
			print STDERR "the format is wrong:",$_;
		}
	}

	close PINF || die;

	print	STDERR "\nread file $pinfile completely,element num=$pk1,$pkj!\n";
	return;
}

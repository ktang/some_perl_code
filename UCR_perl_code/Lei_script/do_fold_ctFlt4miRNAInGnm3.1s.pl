#!/usr/bin/perl
# Copyright (c)  2009-
# Program:			do_fold_ctFlt4miRNAInGnm
# Author:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Program Date:		2009.11.12
# Modifier:			Gaolei <highlei@gmail.com or leigao@ucr.edu>
# Last Modified:	2010.01.05
# Description:	read fasta file, then do unafold, then do ctStrucFilter
#**************************
# Version: 1.1	do blast, blast2histogram, histogram2predict
# Version: 2.0	do unafold for $foldNum sequences each time.
# Version: 2.0s	swinlen	= length of smallRNA
# Version: 3.0s	change the blast to soap and add the tag matching directions into consideration
#**************************
# e-mail:highlei@gmail.com

my $version="3.0s";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #运行开始时间
print STDERR "Now = $Time_Start\n\n";


use Getopt::Std;
getopts("hi:u:d:a:j:m:s:w:c:f:b:n:r:p:l:e:t:0:");
my $flag0		= (defined $opt_0) ? $opt_0 : 0;
my $infile		= $opt_i;
my $fold_set	= (defined $opt_u) ? $opt_u : "-s DAT" ; # --mfold=*,*,5
my $smallRNA	= (defined $opt_d) ? $opt_d : "";
my $multiple4strand	= (defined $opt_a) ? $opt_a : 5;
my $blastORsoap	= (defined $opt_j) ? $opt_j : "blast";

my $mismatch	= (defined $opt_m) ? $opt_m : 4;
my $slideLen	= (defined $opt_s) ? $opt_s : 1;
my $winLen		= (defined $opt_w) ? $opt_w : 21;
my $colInct		= (defined $opt_c) ? $opt_c : 4;	## the pair column
my $maxFree		= (defined $opt_f) ? $opt_f : -35;
my $bulges_size	= (defined $opt_b) ? $opt_b : 2;
my $bulges_num	= (defined $opt_n) ? $opt_n : 2;
my $outRNA		= (defined $opt_r) ? $opt_r : "";
my $path		= (defined $opt_p) ? $opt_p : "./";
my $foldNum		= (defined $opt_l) ? $opt_l : 7000;	# memory: 10,000 ~ 2.5G; 7,000 ~ 2G; 5,000 ~ 1.6G; 1,000 ~ 900m

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

my $colofname	= (defined $opt_a) ? $opt_a : 3;

if ($opt_h){
	usage();
}

my $unafold		= "hybrid-ss-min -s DAT";
#my $blast2hist	= "blast2histogram1.1_m8.pl $b2h_set";
my $hist2predict= "histogram2predict4miRNAInGnm.pl -0 0 $h2p_set";
my ($hsp_1st,$hit_1st,$evalue,$identity,$score,$len,$overlap_len);
my ($minHitsNum,$minRatio);

sub numerically{$a<=>$b};

use FileHandle;
use strict;



my ($i,$j,$k,$num,$ntlen,$k1,$k2,$k3,$m,$n,$file,$line,$line1,$in,$match,$omatch,$a,$b,$end);
my (@bufi,@bufo,@genome,@gnmName,@gnmLen);
my (%ct,%mrna,%gnm,%query);
my $key="";
my $key2="";
my ($ii,$jj,$h2p);

#===========================================================================================================
#====================                  main
#===========================================================================================================
#my $flag0	= 1;
my $yesorno	= "y";
while ($flag0) {
	print STDERR ("\n------------------------------------------------------------\n");
	print STDERR ("\n $0 version $version\n\n");
	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %55s : %-25s","input fasta file name",$infile);#%55s
	printf STDERR ("\n u  %55s : %-25s","the parameter set for unafold",$fold_set);
	printf STDERR ("\n d  %55s : %-25s","input small RNA database (formatdb)",$smallRNA);
	printf STDERR ("\n m  %55s : %-25s","number of mismatch in stem",$mismatch);
	printf STDERR ("\n s  %55s : %-25s","length of slide each time",$slideLen);
	printf STDERR ("\n w  %55s : %-25s","length of window",$winLen);
	printf STDERR ("\n f  %55s : %-25s","allowable max free energy",$maxFree);
	printf STDERR ("\n b  %55s : %-25s","allowable max asymmetric bulge size",$bulges_size);
	printf STDERR ("\n n  %55s : %-25s","allowable max asymmetric bulges num",$bulges_num);
	printf STDERR ("\n r  %55s : %-25s","output the possible mature miRNA sequences",$outRNA);
	printf STDERR ("\n p  %55s : %-25s","input the path of perl scripts",$path);#%55s
	printf STDERR ("\n l  %55s : %-25s","input number of sequcnes for unafold once",$foldNum);#%55s
	printf STDERR ("\n t  %55s : %-25s","the parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n e  %55s : %-25s","the parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
	printf STDERR ("\n x  %55s !","exit the program");

	print STDERR ("\n\n");
	print STDERR "y to accept these or type the letter for one to change!\n";
	$yesorno = <STDIN>;	$yesorno =~s/[\s|\t|\r|\n]+$//g;	$yesorno = lc($yesorno);
	if ($yesorno eq "y") {print STDERR ("\n\n"); $flag0 = 0;}
	elsif($yesorno eq "i") {print STDERR "please input fasta file name:\n"; $infile	= <STDIN>;	$infile	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "u") {print STDERR "please input the parameter set for unafold\n";$fold_set	= <STDIN>;$fold_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "d") {print STDERR "please input small RNA database,should be formatdb:\n"; $smallRNA	= <STDIN>;$smallRNA	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "m") {print STDERR "please input the number of mismatch in stem\n"; $mismatch	= <STDIN>;$mismatch	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "s") {print STDERR "please input length of slide each time\n"; $slideLen	= <STDIN>;$slideLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "w") {print STDERR "please input length of window\n"; $winLen	= <STDIN>;$winLen	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "f") {print STDERR "please input allowable max free energy\n"; $maxFree	= <STDIN>;$maxFree	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "b") {print STDERR "please input allowable max asymmetric bulge size\n"; $bulges_size	= <STDIN>;$bulges_size	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "n") {print STDERR "please input allowable max asymmetric bulge number\n"; $bulges_num	= <STDIN>;$bulges_num	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "r") {print STDERR "please input file name to store the possible mature miRNA sequences\n"; $outRNA	= <STDIN>;$outRNA	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "p") {print STDERR "please input the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl):\n";$path	= <STDIN>;$path	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "l") {print STDERR "please input the number of sequcnes for unafold once:\n";$foldNum	= <STDIN>;$foldNum	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "t") {print STDERR "please input the parameter set for blast2histogram:\n";$b2h_set	= <STDIN>;$b2h_set	=~s/[\s|\t|\r|\n]+$//g;}
	elsif($yesorno eq "e") {print STDERR "please input the parameter set for histogram2predict4miRNAInGnm.pl:\n";$h2p_set	= <STDIN>;$h2p_set	=~s/[\s|\t|\r|\n]+$//g;}

	elsif($yesorno eq "x") {sub_end_program();;exit(0);}
}

#my $minHitsNum	= (defined $opt_e) ? $opt_e : 5;
#my $minRatio	= (defined $opt_a) ? $opt_a : 0.75;
$h2p_set =~ s/\,/ /g;

	print STDERR ("Settings for this run:");
	printf STDERR ("\n i  %55s : %-25s","input fasta file name",$infile);#%55s
	printf STDERR ("\n u  %55s : %-25s","the parameter set for unafold",$fold_set);
	printf STDERR ("\n d  %55s : %-25s","input small RNA database (formatdb)",$smallRNA);
	printf STDERR ("\n m  %55s : %-25s","number of mismatch in stem",$mismatch);
	printf STDERR ("\n s  %55s : %-25s","length of slide each time",$slideLen);
	printf STDERR ("\n w  %55s : %-25s","length of window",$winLen);
	printf STDERR ("\n f  %55s : %-25s","allowable max free energy",$maxFree);
	printf STDERR ("\n b  %55s : %-25s","allowable max asymmetric bulge size",$bulges_size);
	printf STDERR ("\n n  %55s : %-25s","allowable max asymmetric bulges num",$bulges_num);
	printf STDERR ("\n r  %55s : %-25s","output the possible mature miRNA sequences",$outRNA);
	printf STDERR ("\n p  %55s : %-25s","input the path of perl scripts",$path);#%55s
	printf STDERR ("\n l  %55s : %-25s","input number of sequcnes for unafold once",$foldNum);#%55s
	printf STDERR ("\n t  %55s : %-25s","the parameter set for blast2histogram",$b2h_set);
	printf STDERR ("\n e  %55s : %-25s\n","the parameter set for histogram2predict4miRNAInGnm.pl",$h2p_set);
#} else {
#	print STDERR "please input the correct histogram2predict parameter!\n";
#	usage();
#}
if ($b2h_set =~ /([\d|\.|\-]+)\,([\d|\.|\-]+)\,([\d|\.|\-]+)\,([\d|\.|\-]+)\,([\d|\.|\-]+)\,([\d|\.|\-]+)\,([\d|\.|\-]+)/) {
	$hsp_1st	= $1;	$hit_1st	= $2;	$evalue	= $3;	$identity	= $4;
	$score	= $5;	$len	= $6;	$overlap_len	= $7;

	print STDERR "hsp_1st=$hsp_1st,hit_1st=$hit_1st,evalue=$evalue,identity=$identity,score=$score,len=$len,overlap_len=$overlap_len\n\n";
} else {
	print STDERR "please input the correct blast2histogram parameter:$b2h_set!\n";
	usage();
}
$m	= 0;
############################################ read smallRNA files ######################################################
if ($smallRNA ne "") {
	$i	= -1;
	$file = new FileHandle ("$smallRNA") || die("Cannot open the small RNA $smallRNA!\n");
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
	print STDERR "There are $m sequence in small RNA $smallRNA\n";
}

############################################ read genome files ######################################################
$i	= -1;
$file = new FileHandle ("$infile") || die("Cannot open file $infile!\n");
while (<$file>) {
	$_=~s/^[\s|\t]+//g;
	$_=~s/[\s|\t|\r|\n]+$//g;
	if ($_ =~/^>(\S+)/) {
			$i++;
			$genome[$i]	= "";
			$gnmName[$i]	= $_."\n";	#$k1++;
			$query{$1}	= $i;
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
open(HISTF, ">$infile.ctFlt.hist") || die("Can not open file: $infile.ctFlt.hist\n");

$ii	= 0;	$jj	= 0;
$j	= "__do_fold_$infile\_" . $start . ".fa";
for ($k2=0; $k2 < @genome;$k2+=$foldNum) {
	open(TMPF, ">$j") || die("Can not open file: $j\n");
	$k3	= @genome;	$k3 = $k3 > $k2+$foldNum ? $k2+$foldNum : $k3;	print STDERR "\trun to $k3 sequences\n";
	for (my $k5 = $k2; $k5 < $k3 ;$k5++) {
		print TMPF $gnmName[$k5],$genome[$k5],"\n";
	}
	close(TMPF);
	print STDERR "\t  Now = [",sub_format_datetime(localtime(time())),"]\tdo $unafold for $k3 seq.\n";
	system("$unafold $fold_set $j > $j.unafold.out");
	print STDERR "\t  Now = [",sub_format_datetime(localtime(time())),"]\tdo readCTfile for $k3 seq.";
	%ct	= ();
	sub_readCTfile($j.".ct");
	$i	= 0;
	foreach $key2 (sort keys %ct) {
		foreach $key (sort numerically keys %{$ct{$key2}}) {
			if ($key eq "num") {
				next;
			}
			if ($ct{$key2}{$key}{"FE"} > $maxFree) {
				delete($ct{$key2}{$key});
				$ct{$key2}{"num"}--;
				$ii++;
			}
		}
		if ($ct{$key2}{"num"}	== 0) {
			delete($ct{$key2});
	#		next;
		}
	}

	print STDERR "\t\tdelete ",$ii," seq for free energy\n";
#	print STDERR "There are ",$ct{"num"}," sequences remained after check free energy\n\n";

#-----------------------------------------v3.1 blast2histogram, histogram2predict--------------------------------------------
#	system("less $j | gawk \'{if(\$_~/>/){print \$1} else {print \$_}}\' > $j.fna");
	if ($blastORsoap =~/^b/) {
		print STDERR "\t  Now = [",sub_format_datetime(localtime(time())),"]\tdo blast2histogram for ",$k3-$ii," seq.";
#		print STDERR "\t\tblastall -p blastn -F F -v 1000000 -b 1000000 -i $j -d $smallRNA -o $j.blastn -m 8\n";
		system("blastall -p blastn -F F -v 100000000 -b 100000000 -e 1 -i $j -d $smallRNA -o $j.blastn -m 8");
		print STDERR "\t\tblast done!";	#	system("cat $j.blastn > tmp.500.blastn");
		($m,$n,$h2p) =	sub_blast2histogram("$j.blastn");	print STDERR "\tn=$n\n";
	} else {
		print STDERR "\t  Now = [",sub_format_datetime(localtime(time())),"]\tdo soap2histogram for ",$k3-$ii," seq.";
		system("2bwt-builder $j 2> $j.soap.err");
		$m	= `cat $j.soap.err`;
		if ($m ne "") {
			print STDERR "Wrong when do 2bwt-builder!\n";
			last;
		}
#		print STDERR "\t\tsoap -D $j.index -M 0 -r 2 -v 0 -a $smallRNA -o $j.soap 2> $j.soap.err\n";
		system("soap -D $j.index -M 0 -r 2 -v 0 -a $smallRNA -o $j.soap 2> $j.soap.err");
		print STDERR "\t\tsoap done!";
#		system("cat $j.soap > tmp.500.soap");
	#	$m	= $genome[$k2];	$m=~s/[\s|\t|\r|\n]+//g;
		($m,$n,$h2p) =	sub_soap2histogram("$j.soap");	print STDERR "\tn=$n\n";#h2p=\n$h2p";
	}
	if ($n	== 0) {
#		print STDERR "-------here n=$n\n";
		next;
	}
	print STDERR "\t  Now = [",sub_format_datetime(localtime(time())),"]\tdo strucCheck for ",$k3-$ii," seq.\n";
#-----------------------------------------v3.1 blast2histogram, histogram2predict--------------------------------------------
	($i,$a,$b,$file)	= sub_strucCheck();			# i=num, a=ctflt.txt, b=stem, file=ctfile.
	if ($i	== 0) {
		next;
	}
	$jj	+= $i;
	print	CTFLT	$a;		#,"\n";
	print	STEMF	$b;		#,"\n";
	print	CTF		$file;	#,"\n";
	print	HISTF	$h2p;	#,"\n";

#	$i	= <STDIN>;
}
close(CTFLT) || die("Can not close file: $infile.ctFlt.txt\n");
close(STEMF) || die("Can not close file: $infile.ctFlt.stem\n");
close(CTF) || die("Can not close file: $infile.ctFlt.ct\n");
close(HISTF) || die("Can not close file: $infile.ctFlt.hist\n");

print STDERR "perl $path/$hist2predict\n";
system("perl $path/$hist2predict -i $infile.ctFlt.hist -c $infile.ctFlt.txt > $infile.ctFlt.h2p 2>> $j.err.h2p");
system("cat $infile.ctFlt.h2p >> $j.out.h2p");
print STDERR "\nThere are ",$ii," sequences which are deleted due to free energy > $maxFree.\n";
print STDERR "There are ",$jj," sequences remained after structure check.\n\n";
system("rm $j*");



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
	print "			eg: -u --mfold=*,*,5; default: -s DAT\n";
	print "-d	input the small RNA database,should be formatdb\n";
	print "			eg: smallRNA.fa\n";
	print "-m	input number of mismatch in stem\n";
	print "			eg: -m 5; default: 4\n";
	print "-s	input length of slide each time\n";
	print "			eg: -s 2; default: 1\n";
	print "-w	please input length of window\n";
	print "			eg: -w 22; default: 21\n";
	print "-f	please input allowable max free energy\n";
	print "			eg: -f -40; default: -35\n";
	print "-b	please input allowable max asymmetric bulge size\n";
	print "			eg: -b 1; default: 2\n";
	print "-n	please input allowable max asymmetric bulges number\n";
	print "			eg: -n 0; default: 2\n";
	print "-r	output the possible mature miRNA sequences\n";
	print "			eg: -o possible_miRNA.fa; default: donot output\n";
	print "-p	input the path of perl scripts (qsub_noFasta.pl & ctStructureFilter.pl)\n";
	print "			eg: -p ../; default: ./\n";
	print "-l	input the number of sequcnes for unafold once\n";
	print "			eg: -l 10,000; default: 2,000\n";
	print "-t	the parameter set for blast2histogram.pl\n";
	print "			eg: -e \"\"; default: $b2h_set\n";
	print "-e	the parameter set for histogram2predict4miRNAInGnm.pl\n";
	print "			eg: -t \"\"; default: $h2p_set\n";
	print "-h	display this lines\n";
	print "\nExample:\n";
	print "$0 -i RNAfold.ct -m 4\n";
	print "$0 -i RNAfold.ct -m 5 -s 1 -w 21\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 2\n";
	print "$0 -i RNAfold.ct -m 4 -s 1 -w 21 -f -35 -b 2 -n 2 -r possible.miRNA.fa\n";
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
######################                  sub_readCTfile
############################################################################################################

sub sub_readCTfile
{
	my ($rcname) = @_;
	my ($rcj,$rci,$rck,$rcm,$rcn,$rcfile);
	my (@rcbuf);

	$rcfile = new FileHandle ("$rcname") || die("Cannot open ct file: $rcname!\n");
	while (<$rcfile>) {
		$_=~s/^[\s|\t]+//g;
		$_=~s/[\s|\t|\r|\n]+$//g;
		splice(@rcbuf,0);
		if ($_ =~/\=/) {
			$rcj			= $_;
			@rcbuf			= split(/\t+/,$rcj);
#			$rck = "";
#			for ($rcm = $colofname; $rcm <= @rcbuf;$rcm++) {
#				$rck	.= $rcbuf[$rcm-1];
#			}
			$rck	=	$rcbuf[$colofname - 1];	
			if (!exists($ct{$rcbuf[0]}{"num"})) {
				$ct{$rck}{"num"}		= 1;	$rcn	= 1;
			} else {
				$ct{$rck}{"num"}++;	$rcn	= $ct{$rck}{"num"};	#print STDERR "n=$n,num=",$ct{$k}{"num"},"\n";
			}
			$ct{$rck}{$rcn}{"name"}		= $rcbuf[2];		## name of sequence
			$ct{$rck}{$rcn}{"len"}		= $rcbuf[0];
			$ct{$rck}{$rcn}{"title"}		= $rcj;		#	print STDERR $ct{$rcn}{"title"},"=title in read ct file\n";
			$rcj			=~ /\=\s?(\S+)/;
			$ct{$rck}{$rcn}{"FE"}		= $1;
			$rci			= 0;
		} else {
			$rcj	= $_;
			@rcbuf    = split(/\s+/,$rcj);
			for ($rcm = 0; $rcm < @rcbuf;$rcm++) {
				$ct{$rck}{$rcn}{$rci}{$rcm}	= $rcbuf[$rcm];
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
	my ($sk,$sj,$skey,$skey2,$sk1,$sk2,$sm,$sn,$sline,$sline1,$sflag,$sflag0,$opRNA,$opRNA0,$sct);	#,$ovlp_s,$ovlp_e
	my (@sbuf,@stmp,@stmp1);
	my %smrna=();
	my ($ovlp_s,$ovlp_e,$swinLen);

	$sline	= "";	$opRNA0	= "";	$sj	= 0;	$sct	= "";
foreach $skey2 (sort keys %ct) {

	foreach $skey (sort numerically keys %{$ct{$skey2}}) {
		if ($skey eq "num") {
			next;
		}
		$ct{$skey2}{$skey}{"title"}	=~ /(\d+)\,(\d+)\s*$/;
		$ovlp_s	= $1;	$ovlp_e	= $2;	#print STDERR $ct{$skey}{"title"},"\ts=$ovlp_s\te=$ovlp_e\n";
		$swinLen	= $ovlp_e - $ovlp_s + 1;# > $winLen ? $winLen : $ovlp_e - $ovlp_s + 1;	#print STDERR "swinLen=$swinLen\n";
		splice(@sbuf,0);				## the pair column
		foreach $sk (sort numerically keys %{$ct{$skey2}{$skey}}) {
			if ($sk eq "title" || $sk eq "FE" || $sk eq "name" || $sk eq "len") {
				next;
			}
			$sbuf[$sk]	= $ct{$skey2}{$skey}{$sk}{$colInct};
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
			$opRNA	.= 	">".$ct{$skey2}{$skey}{"name"}."_" . ($sk1+1). ".." . ($sk1+$swinLen) . "\n";
			for ($sm = $sk1;$sm < $sk1+$swinLen ;$sm++) {
				$opRNA	.= $ct{$skey2}{$skey}{$sm}{1};
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
			$sline	.= "\n".$ct{$skey2}{$skey}{"name"}."\t" . $sline1;
			$sct	.= $ct{$skey2}{$skey}{"title"}."\n";	#print STDERR $ct{$skey}{"title"},"=title,ct=$sct----\n";
			foreach $sk (sort numerically keys %{$ct{$skey2}{$skey}}) {
				if ($sk eq "title" || $sk eq "FE" || $sk eq "name" || $sk eq "len") {
					next;
				}
				foreach $sm (sort numerically keys %{$ct{$skey2}{$skey}{$sk}}) {
					$sct	.= $ct{$skey2}{$skey}{$sk}{$sm}."\t";
				}
		#		$sct	.= "\n";
				$sct	= substr($sct,0,length($sct)-1)."\n";
			}
		#	print STDERR "ct=$sct\n";
		}
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

############################################################################################################
######################                  sub_blast2histogram
############################################################################################################

sub sub_blast2histogram 
{
	my ($binfile) = @_;
	my ($bm,$bn,$bk1,$bk2,$bk3,$bk4,$ba,$bb,$bii,$bkk,$bjj,$bflag,$bfile,$bop,$bstrand,$bnum);
	my ($qlen,$dblen,$qname, $dbname, $qidentity, $alnLen, $qstart, $qend, $dbstart, $dbend, $qevalue, $qbits);
	my (%entry,%bmorp,%bbuf);
	my (@bbufm,@bbufp);

#	if ($outfile ne "") {
#		open(OUT, ">$outfile") || die("Can not open file: $outfile\n");
#	}
	$bk1	= 0;	$bflag	= 0;	$ba = 0;	$bm = 0;	$bn	= 0;	$bop = "";	$bnum = 0;
	$bfile = new FileHandle ("$binfile") || die("Cannot open file: $binfile");

	while(<$bfile>)
	{
		if ($_ =~/^(\S+)\s+.+\s+(\S+)\s+([\d|\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([\d|\.]+)/ ||
			$_ =~/^(\S+)\s+(\S+)\s+([\d|\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([\d|\.]+)/) {
			next if $9 >= $evalue;	$bb	= $_;
			if (!exists($entry{$1})) {
				%entry=();
		#		print STDERR "\n\tB $qname:\t+:",$bmorp{"+"},"\t\*$multiple4strand:",$multiple4strand*$bmorp{"+"},"\t-:",$bmorp{"-"},"";

				if ($bmorp{"+"} >= $multiple4strand*$bmorp{"-"}) {
					$bstrand	= "+";
				} elsif ($bmorp{"-"} >= $multiple4strand*$bmorp{"+"}) {
					$bstrand	= "-";
				} else {
					$bmorp{"-"} = 0;	$bmorp{"+"} = 0;	$bstrand	= "+";
				}
				if ($bmorp{$bstrand} != 0) {
					$bop .= ">".$qname."\n";
					for ($bk2 = 0;$bk2 < $qlen ;$bk2++) {
						$bop .= $bbuf{$bstrand}->[$bk2]." ";
					}
					$bop .= "\n";	$bnum++;
				} else {
	#				delete($entry{$qname});
	#				print STDERR "a=$ba---\n";
				}

				$entry{$1}{"num"}	= 1;#	$bflag	= 0;
				$qlen	= length($genome[$query{$1}] );
				$dblen	= length($gnm{$2});	%bbuf=();	#	splice(@bbuf,0);	splice(@bbufm,0);	splice(@bbufp,0);
				for ($bk2 = 0; $bk2 < $qlen ;$bk2++) {
					$bbuf{"+"}->[$bk2]	= 0;	$bbuf{"-"}->[$bk2]	= 0;#	$bbufm[$bk2]	= 0;	$bbufp[$bk2]	= 0;
				}
				$bkk	= 0;	$bk1	= 0;	$bmorp{"-"} = 0;	$bmorp{"+"} = 0;	$bii++;	#$ba	= 0;	
			} else {
				$dblen	= length($gnm{$2});
				if ($hit_1st != -1 && $bkk > $hit_1st) { 
	#				print STDERR "kk=$bkk---\n";
					next;
				}
			}
			$qname=$1; $dbname=$2; $qidentity=$3; $alnLen=$4; $qevalue=$9; $qbits=$10;
			$dbstart	= $7 > $8 ? $8 : $7;	$dbend		= $7 > $8 ? $7 : $8;
			$qstart		= $5 > $6 ? $6 : $5;	$qend		= $5 > $6 ? $5 : $6;
			$bstrand	= $7 > $8 ? "-" : "+";
	#		$qstart=$5; $qend=$6; $dbstart=$7; $dbend=$8; 
			
			$dbname	=~ /\_(\d+)x?$/;
			$bk3	= $1;
			if ($qevalue	== 0 && $qbits < 50) {
				$bkk--;	#print STDERR "qevalue=$qevalue---\n";
				next;
			}
			if (exists($entry{$qname}{$dbname})) {
				$bjj++;
				$entry{$qname}{$dbname}++;
				if ($hsp_1st == 0) {
				} elsif ($entry{$qname}{$dbname} > $hsp_1st){ #$hsp_1st == 1 && $bjj > 1) {
	#				print STDERR "jj=$bjj---\n";
					next;
				}
			} else {
				$bjj	= 1;
				$entry{$qname}{$dbname}	= 1;
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
			if ($bkk == 0) {
				$bn++;	$bkk++;
			} else {
				$bkk++;
			}
			$bmorp{$bstrand}	+= $bk3;
			for ($bk2 = $qstart -1; $bk2 < $qend ;$bk2++) {
				$bbuf{$bstrand}->[$bk2]	+= $bk3;	$ba	+= $bk3;
			}
	#		if ($bjj == 1 && $outfile ne "") {
	#			print OUT $qname,"\t",$dbname,"\n";
	#			if ($database ne "") {
	#				if (exists($gnm{$dbname})) {
	#					print OUT $gnm{$dbname};
	#				} else {
	#					print STDERR "no seq:",$dbname,"\t",$qname,"\n";
	#				}
	#			}
	#		}
		$bm++;	#	print "$dbname\t$qname\t$bstrand\n";
		} else {
			print STDERR "the format is wrong at $bm:",$_;
		}
	}
#	print STDERR "\n\tB $qname:\t+:",$bmorp{"+"},"\t\*$multiple4strand:",$multiple4strand*$bmorp{"+"},"\t-:",$bmorp{"-"},"";

	if ($bmorp{"+"} >= $multiple4strand*$bmorp{"-"}) {
		$bstrand	= "+";
	} elsif ($bmorp{"-"} >= $multiple4strand*$bmorp{"+"}) {
		$bstrand	= "-";
	} else {
		$bmorp{"-"} = 0;	$bmorp{"+"} = 0;	$bstrand	= "+";
	}
	if ($bmorp{$bstrand} != 0) {
		$bop .= ">".$qname."\n";
		for ($bk2 = 0;$bk2 < $qlen ;$bk2++) {
			$bop .= $bbuf{$bstrand}->[$bk2]." ";
		}
		$bop .= "\n";	$bnum++;
	}
#	if ($ba != 0) {
#		$bop .= ">".$qname."\n";
#		for ($bk2 = 0;$bk2 < $qlen ;$bk2++) {
#			$bop .= $bbuf[$bk2]." ";
#		}
#		$bop .= "\n";
#	}

	close $bfile || die;
	return $bn,$bnum,$bop;
}

############################################################################################################
######################                  sub_soap2histogram
############################################################################################################

sub sub_soap2histogram 
{
	my ($binfile) = @_;
	my ($bm,$bn,$bk1,$bk2,$bk3,$bk4,$ba,$bb,$bii,$bkk,$bjj,$bflag,$bfile,$bop);
	my ($qlen,$dblen,$qname, $dbname, $qidentity, $alnLen, $qstart, $qend, $dbstart, $dbend, $qevalue, $qbits);
	my (%entry);
	my (@bbuf);

#	if ($outfile ne "") {
#		open(OUT, ">$outfile") || die("Can not open file: $outfile\n");
#	}
	$bk1	= 0;	$bflag	= 0;	$ba = 0;	$bm = 0;	$bn	= 0;	$bop = "";
	$bfile = new FileHandle ("$binfile") || die("Cannot open file: $binfile");
	while (<$bfile>) {
		if ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+|-]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+/) {
			if ($10 != 0) {
	#			print STDERR "\tMM!=0,",$_;
				next;
			}
			if (exists($entry{$8})) {
			} else {
				$entry{$8}{"+"}	= 0;
				$entry{$8}{"-"}	= 0;	#print STDERR "$7,";
			}
	#		$entry{$8}{$7}++;
			$entry{$8}{$1}->[0]	= $7;			# +/-
			$entry{$8}{$1}->[1]	= $9;			# start
			$entry{$8}{$1}->[2]	= $9+$6-1;		# end
#			$entry{$8}{$1}->[3]	= $6;			# Len
			$bk2	= $1;	$bk1	= $8;	$bk4	= $7;
			$bk2	=~  /\_(\d+)x?$/;
		#	$bk3	= $1;
			$entry{$bk1}{$bk4}	+=$1;
		} else {
			print STDERR "In soap2hist, format is wrong:",$_;
		}
	}
	close $bfile || die;
	foreach $bk1 (sort keys %entry) {	
#		print STDERR "\n\tS $bk1:\t+:",$entry{$bk1}{"+"},"\t\*$multiple4strand:",$multiple4strand*$entry{$bk1}{"+"},"\t-:",$entry{$bk1}{"-"},"";
		if ($entry{$bk1}{"+"} >= $multiple4strand*$entry{$bk1}{"-"}) {
			$bflag	= "+";
		} elsif ($entry{$bk1}{"-"} >= $multiple4strand*$entry{$bk1}{"+"}) {
			$bflag	= "-";
		} else {
			next;
		}
		$qlen	= length($genome[$query{$bk1}] );
		splice(@bbuf,0);
		for ($bk2 = 0; $bk2 < $qlen ;$bk2++) {
			$bbuf[$bk2]	= 0;
		}
		$ba = 0;	
		foreach $bk2 (sort keys %{$entry{$bk1}}) {
			if ($bk2 eq "+" || $bk2 eq "-") {
				next;
			}
			if ($entry{$bk1}{$bk2}->[0] ne $bflag) {
				next;
			}
			$bk2	=~  /\_(\d+)x?$/;
			$bk3	= $1;
			for ($bk4 = $entry{$bk1}{$bk2}->[1]-1;$bk4 < $entry{$bk1}{$bk2}->[2] ;$bk4++) {
				$bbuf[$bk4]	+= $bk3;	$ba	+= $bk3;
			}
		}
		if ($ba != 0) {
			$bop .= ">".$bk1."\n";
			for ($bk2 = 0;$bk2 < $qlen ;$bk2++) {
				$bop .= $bbuf[$bk2]." ";
			}
			$bop .= "\n";	$bm++;
		}
	}
	%entry=();
#	close $bfile || die;
	return $bn,$bm,$bop;
}
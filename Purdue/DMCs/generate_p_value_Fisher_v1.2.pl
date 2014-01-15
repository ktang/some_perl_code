#! /usr/bin/perl -w
#		   WT	    mut
#		  word2   ~word2
#mC	 	 word1    n11      n12 | n1p
#No.T	~word1    n21      n22 | n2p
#          --------------
#		  np1      np2   npp
use strict;
use Text::NSP::Measures::2D::Fisher2::twotailed;

#v1.0 use Dr. Liu's output as input to prove the correctness
#of the code.

#v1.1 use Col-0 high WT

#v1.2 use C24 input

#my ($npp, $n1p,    $np1,      $n11);
#   total  mC_sum   WT_depth  WT_mC;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $wt_f1 = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/cal_methy/onetype_chr/JKZ19_C24_Luc_149_s2_meth_db_input_onetype_chr.txt";
my $wt_f2 = "/Volumes/Macintosh_HD_2/BS_data_HiSeq/C24_background/cal_methy/onetype_chr/JKZ20_C24_Luc_150_s3_meth_db_input_onetype_chr.txt";

my $usage = "$0 <mut_file1> <mut_file2> <output>";

die $usage unless (@ARGV == 3);
#my $input = $ARGV[0];

my $mut_f1 = $ARGV[0];
my $mut_f2 = $ARGV[1];

my $output = $ARGV[2];

if($debug){
	$wt_f1  = "/Users/tang58/try/debug/fisher/JKZ19_C24_Luc_149_s2_meth_db_input_onetype_chr_1000line.txt";
	$wt_f2  = "/Users/tang58/try/debug/fisher/JKZ20_C24_Luc_150_s3_meth_db_input_onetype_chr_1000line.txt";
	$mut_f1 = "/Users/tang58/try/debug/fisher/JKZ23_ros3_1_153_s7_meth_db_input_onetype_chr_1000line.txt";
	$mut_f2 = "/Users/tang58/try/debug/fisher/JKZ24_ros3_1_154_s8_meth_db_input_onetype_chr_1000line.txt";
}

die "wrong input" unless (-e $wt_f1 and -e $wt_f2 and -e $mut_f1 and -e $mut_f2);
die "output exists" unless (! (-e $output));

print STDERR "WT:\n$wt_f1\n$wt_f2\n\n";
print STDERR "mut:\n$mut_f1\n$mut_f2\n\n";



#open (IN, $input) or die "cannot open $input: $!";

open (WT1, $wt_f1) or die "cannot open $wt_f1 : $!";
open (WT2, $wt_f2) or die "cannot open $wt_f2: $!";
open (MUT1, $mut_f1) or die "cannot open $mut_f1: $!";
open (MUT2, $mut_f2) or die "cannot open $mut_f2 : $!";

open (OUT, ">$output") or die "cannot open $output: $!";

my ($l_wt1, $l_wt2, $l_mut1, $l_mut2);


#chr1	1	0		0	0	+		CHH		0
# 0     1   2		3	4	5		6		7
#chr	pos	depth	mC	per	strand	type 	isMet
while ($l_wt1 = <WT1>, $l_wt2 = <WT2>, $l_mut1 = <MUT1>, $l_mut2 = <MUT2>){
 	chomp $l_wt1;
 	chomp $l_wt2;
 	chomp $l_mut1;
 	chomp $l_mut2;
 	
 	my @a_wt1 = split "\t", $l_wt1;
 	my @a_wt2 = split "\t", $l_wt2;
 	my @a_mut1 = split "\t", $l_mut1;
 	my @a_mut2 = split "\t", $l_mut2;
 	
 	die "$l_wt1, $l_wt2, $l_mut1, $l_mut2" unless ( $a_wt2[0] eq $a_wt1[0] and $a_mut1[0] eq $a_wt1[0] and $a_mut2[0] eq $a_wt1[0] and 
 	 $a_wt2[1] == $a_wt1[1] and $a_mut1[1] == $a_wt1[1] and $a_mut2[1] == $a_wt1[1]);
 
 
    my ($npp, $n1p,    $np1,      $n11);
#    chomp;
#    my @a = split "\t";
#    my $wt_depth = $a[4];
#    my $mut_depth = $a[7] + $a[10];
	my $wt_depth  = $a_wt1[2]  + $a_wt2[2];
	my $mut_depth = $a_mut1[2] + $a_mut2[2];
	
	
	if($wt_depth >= 4 and $mut_depth >= 4) {
#    if($wt_depth >= 4 and $mut_depth >= 4 and $wt_depth <= 100  and $mut_depth <= 100 ){
#		my $wt_mC = round($wt_depth * $a[5] );
#		my $mut_mC = round($a[7] * $a[8]) + round($a[10] * $a[11]);

		my $wt_mC  = $a_wt1[3] + $a_wt2[3];
		my $mut_mC = $a_mut1[3] + $a_mut2[3];
		
		$npp = $wt_depth + $mut_depth;
		$n1p = $wt_mC + $mut_mC;
		$np1 = $wt_depth;
		$n11 = $wt_mC;
		my  $p_value = calculateStatistic( n11=>$n11,
										   n1p=>$n1p,
						    			   np1=>$np1,
						    			   npp=>$npp);
						    			   
		my $errorCode;
		if( ($errorCode = getErrorCode())){
			print STDERR $errorCode." - ".getErrorMessage();
			die "$l_wt1, $l_wt2, $l_mut1, $l_mut2";
  		}
  		else{
#    print getStatisticName."value for bigram is ".$twotailed_value;
			print OUT join("\t", ($a_wt1[0], $a_wt1[1], $a_wt1[5], $a_wt1[6], @a_wt1[2..4], $a_wt1[7], @a_wt2[2..4], $a_wt2[7], 
			@a_mut1[2..4], $a_mut1[7], @a_mut2[2..4], $a_mut2[7], $p_value)), "\n";

  		}
    }
}

exit;

sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}

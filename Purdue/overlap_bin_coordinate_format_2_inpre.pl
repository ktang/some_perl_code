#!/usr/bin/perl -w
# 1:1-100

use strict;
use File::Spec;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/overlap_bin_coordinate_format.pl";
die unless (-e $script);

my $debug = 0;
if ($debug){
	print STDERR "debug: $debug\n";
}

print STDERR "input must have a head\n\n";
my $usage = "$0 \n <indir1> <pre1> <indir2> <pre2> \n\n";
die $usage unless ( @ARGV == 4 );

my $indir1 = shift or die;
my $inpre1 = shift or die;

my $indir2 = shift or die;
my $inpre2 = shift or die;

die unless (-d $indir1);
die unless (-d $indir2);

my @labels = ( "CG_hyper", "CHG_hyper", "CHH_hyper","CG_hypo", "CHG_hypo",  "CHH_hypo" );

#die unless (-e $input1);
#die unless (-e $input2);

foreach my $label (@labels){
	opendir(DIR1, $indir1) or die;
	opendir(DIR2, $indir2) or die;
	
	my @inputs1 = grep /$inpre1/, grep /$label/, readdir DIR1;
	my @inputs2 = grep /$inpre2/, grep /$label/, readdir DIR2;
	
	
#	print STDERR "1\n";
#	print STDERR join("\n", @inputs1), "\n";
#	print STDERR "2\n";
#	print STDERR join("\n", @inputs2), "\n\n";
	
	die unless (  @inputs1 == 1); 
	die unless (  @inputs2 == 1); 
	
	closedir DIR1;
	closedir DIR2;
	
	my $in1 = File::Spec->catfile ($indir1, $inputs1[0]);
	my $in2 = File::Spec->catfile ($indir2, $inputs2[0]);
	
	my $cmd = " perl $script $in1 $in2";
	if($debug){
		print STDERR $cmd, "\n\n";
	}elsif(!$debug){
		`$cmd`;
	}
}
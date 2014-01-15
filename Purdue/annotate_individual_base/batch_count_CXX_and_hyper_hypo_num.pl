#! /usr/bin/perl -w

#/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Oct_1/P001
#use File::Spec;
#$x= File::Spec->catfile('a', 'b', 'c');

use strict;
use File::Spec;

my $script = "/Volumes/Macintosh_HD_2/DMC_base_analysis/src/count_CXX_and_hyper_hypo_num.pl";
die "script" unless (-e $script);

my $usage = "$0 <indir>";

die $usage unless (@ARGV == 1);

my $indir = $ARGV[0];

opendir(DIR, $indir) or die "cannot open $indir: $!";

my @files = grep /txt/, readdir DIR;

for my $file (@files){
	my $input = File::Spec->catfile($indir, $file);
	
	my $cmd = "perl $script $input";
	`$cmd`;
}
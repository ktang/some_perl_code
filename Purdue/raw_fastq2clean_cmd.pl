#!/usr/bin/perl -w

use strict;

my $debug = 1;

my $trim_script = "/Users/tang58/scripts_all/perl_code/Purdue/raw_fastq_to_clean_fastq_v1.0.pl";

my $adapter_P1 = "AGATCGGAAGAGCACA";

my $adapter_P2 = "AGATCGGAAGAGCGTC";

my $phred = 30;
my $minLen = 24;
my $dir = ".";

my $usage ="$0 <in1> <in2> <pre>";

die $usage unless (@ARGV == 3);

my ($in1, $in2, $pre) = @ARGV[0..2];

my $log = $pre."_trim_log.txt";

my $cmd = "time perl $trim_script --adaptor_P1 $adapter_P1 --adaptor_P2 $adapter_P2 --input_P1 $in1 --input_P2 $in2 --output_pre $pre --dir $dir --minLen $minLen --phred $phred > $log";

print STDERR $cmd, "\n\n";

if(!$debug){
	`$cmd`;
}
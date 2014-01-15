#!/usr/bin/perl -w

#v0.1
# v0.0 is use for 1 tfile and 1 cfile

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $genome_size = "1.2e8";
my $keep_dup = "all";

my $usage = "\n$0 \n\n <outdir> <outpre> <tfile_num> <treatment_ChIP> [...] <cfile_num> <control_Input> [...]\n\n";
die $usage unless(@ARGV >= 6);

#my $tfile     = shift or die;
#my $cfile     = shift or die;
my $outdir    = shift or die;
my $outpre    = shift or die;

my $tfile_num = shift or die;
my @tfiles;
for my $i(1..$tfile_num){
	my $tmp = shift or die;
	push @tfiles, $tmp;
}
my $cfile_num = shift or die;

my @cfiles;
for my $i(1..$cfile_num){
	my $tmp = shift or die;
	push @cfiles, $tmp;
}

my $log_file = File::Spec->catfile($outdir, $outpre . "_MACS2_log.txt");
die "log file exist\n\n" if(-e $log_file);

my $cmd = "time macs2 callpeak -t @tfiles -c @cfiles -g $genome_size --keep-dup $keep_dup --outdir $outdir -n $outpre 2>> $log_file";

print STDERR $cmd, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd, "\n";
	close OUT;
	`$cmd`;
}
exit;
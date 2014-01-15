#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/ChIP_Seq/simplify_peak_xls_files.pl";
die "script" unless (-e $script);

my $usage = "\n$0 \n <indir> <outdir>\n\n";

die $usage unless (@ARGV == 2);
my $indir = shift or die;
die "indir is not a dir" unless (-d $indir);

my $outdir = shift or die;
die "outdir is not a dir" unless(-d $outdir);

opendir(DIR, $indir) or die "cannot open indir";

my @files = grep /\.xls$/, readdir DIR;

foreach my $file (@files){
	if ($file =~ /(\S+)\.xls$/){
		#my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		my $output = File::Spec->catfile($outdir, $1. ".txt");
		my $cmd = " $script $input $output";
		print STDERR $cmd, "\n\n";
		unless($debug){
			`$cmd`;
		}
	}
}
#my $usage = "\n$0 \n\n<indir> <inpre> <outdir> <outpre> <phred_33_64>\n\n";
exit;
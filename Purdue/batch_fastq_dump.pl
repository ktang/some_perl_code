#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
my $debug = 0;
my $usage = "$0 <do>";
die $usage unless(@ARGV == 1 and $ARGV[0] eq "do");

my $dump = "/Users/tang58/Software/NCBI/sratoolkit.2.1.10-mac64/bin/fastq-dump";

my $dir = ".";
opendir(DIR, $dir) or die;
my @files = grep /\.sra$/, readdir DIR;

print STDERR join("\n", @files), "\n\n";

foreach my $input(@files){
	print STDERR "extract $input...\n";
	my $cmd = "$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;

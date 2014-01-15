#!/usr/bin/perl -w

#v0.1
# output gz files

# XXX.fq.gz

use strict;
use File::Spec;


my $debug = 0;
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);
#die $usage unless(@ARGV == 1 and $ARGV[0] eq "do");

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir );
die unless (-d $outdir);

my $dump = "/Users/tang58/Software/NCBI/sratoolkit.2.1.10-mac64/bin/fastq-dump";

my $dir = $indir; #".";
opendir(DIR, $dir) or die;
my @files = grep /\.sra$/, readdir DIR;
close DIR;

print STDERR join("\n", @files), "\n\n";

foreach my $file(@files){
	
	print STDERR "extract $file...\n";
	
	my $input = File::Spec->catfile( $indir, $file);
	my $output = "";
	if ($file =~ /(\S+)\.sra$/){
		$output = File::Spec->catfile( $outdir, $1. ".fq.gz");
	}
	die $input unless (-e $input);
	die $output if( -e $output );
	
	my $cmd = "$dump --split-3 -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input --outdir $outdir --gzip";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;

#!/usr/bin/perl -w

use strict;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}
my %phred_label = ( 33=> "--phred33", 64=>"--phred64");

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/bowtie2/Xiaohong_reseq_Dec24/step1_map_using_bowtie2_PE_Dec24_no_direct_change.pl";
die "script" unless (-e $script);

my $usage = "\n$0 \n <indir> <outdir> <phred_33_64>\n\n";

die $usage unless (@ARGV == 3);
my $indir = shift or die;
die "indir is not a dir" unless (-d $indir);

my $outdir = shift or die;
die "outdir is not a dir" unless(-d $outdir);

my $phred_num = shift or die;
die unless (defined $phred_label{$phred_num});

opendir(DIR, $indir) or die "cannot open indir";

my @files = grep /R1\.fq\.gz$/, readdir DIR;

foreach my $file (@files){
	if ($file =~ /(\S+)_R1\.fq.gz$/){
		my $pre = $1;
		my $cmd = "time $script $indir $pre $outdir $pre $phred_num";
		print STDERR $cmd, "\n\n";
		unless($debug){
			`$cmd`;
		}
	}
}
#my $usage = "\n$0 \n\n<indir> <inpre> <outdir> <outpre> <phred_33_64>\n\n";
exit;
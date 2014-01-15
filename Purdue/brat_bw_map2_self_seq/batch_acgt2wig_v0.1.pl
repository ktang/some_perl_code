#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/brat_bw_map2_self_seq/acgt2wig_v0.1.pl";#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;

die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir) or die "cannot open $indir: $!";
my @files = grep /_rev\.txt$/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	if ( $file =~ /(\S+)_rev\.txt$/) {
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $cmd = "$script $indir $pre $outdir $pre $pre";
		print STDERR $cmd, "\n\n";
		if (!$debug) {
			`$cmd`;
		}
#my $usage = "$0 \n <acgt_dir> <acgt_pre> <outdir> <sample_label_in_wig> <wig_pre_name> \n\n";
		
	}
	
}

exit;



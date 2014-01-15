#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/src/513_downstream/single_sample_count_feature_number_v0.1.pl";#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <output>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
#my $outdir = shift or die;
my $output = shift or die;
die unless (-d $indir);
#die unless (-d $outdir);

opendir(DIR, $indir) or die "cannot open $indir: $!";

my @files = grep /\.txt$/, readdir DIR;

closedir DIR;

foreach my $file(@files){
	my $input = File::Spec->catfile($indir, $file);
	die unless (-e $input);
	
	my $cmd = "perl $script $input >> $output";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;
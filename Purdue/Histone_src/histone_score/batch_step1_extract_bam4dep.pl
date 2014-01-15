#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/step1_extract_bam4dep.pl";
#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);
#<input_list_DMR> <bam_dir> <outdir> <output_pre>
if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <DMR_indir> <bam_dir> <outdir>\n\n";
die $usage unless(@ARGV == 3);

my $indir_DMR = shift or die;
my $bam_dir = shift or die;
my $outdir = shift or die;

#die unless (-d $indir);
die unless (-d $indir_DMR);
die unless (-d $bam_dir);
die unless (-d $outdir);

opendir(DIR, $indir_DMR) or die "cannot open $indir_DMR: $!";
my @DMR_lists = grep /\.txt$/ , readdir DIR;
closedir DIR;

if ( $debug ) {
	print STDERR join("\n", @DMR_lists), "\n\n";
}

foreach my $file(@DMR_lists){
	if ( $file =~ /(\S+)\.txt$/) {
		my $pre = $1;
		my $input = File::Spec->catfile( $indir_DMR, $file );
		die unless (-e $input);
		my $cmd = "perl $script $input $bam_dir $outdir $pre";
		print STDERR $cmd, "\n\n";
		if ( !$debug) {
			`$cmd`;#code
		}
	}
}
# <input_list_DMR> <bam_dir> <outdir> <output_pre>
exit;

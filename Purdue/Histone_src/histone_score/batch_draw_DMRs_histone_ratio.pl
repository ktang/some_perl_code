#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/draw_DMRs_histone_ratio.r";
#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);
#<input_list_DMR> <bam_dir> <outdir> <output_pre>
if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir_DMR = shift or die;
#my $bam_dir = shift or die;
my $outdir = shift or die;

#die unless (-d $indir);
die unless (-d $indir_DMR);
#die unless (-d $bam_dir);
die unless (-d $outdir);

opendir(DIR, $indir_DMR) or die "cannot open $indir_DMR: $!";
my @DMR_lists = grep /\.txt$/ , readdir DIR;
closedir DIR;

if ( $debug ) {
	print STDERR join("\n", @DMR_lists), "\n\n";
}

my $i = 0;
foreach my $file(@DMR_lists){
	if ( $file =~ /(\S+)\.txt$/) {
		$i++;
		my $pre = $1;
		my $input = File::Spec->catfile( $indir_DMR, $file );
		die unless (-e $input);
		
		my $cmd = "R --slave --vanilla --args $input $pre $outdir < $script";
#		my $cmd = "perl $script $input $bam_dir $outdir $pre";
		
		if ( !$debug) {
			print STDERR $i, "\t";
			`$cmd`;#code
		}else{
			print STDERR $cmd, "\n\n";
		}
	}
}
# <input_list_DMR> <bam_dir> <outdir> <output_pre>
exit;

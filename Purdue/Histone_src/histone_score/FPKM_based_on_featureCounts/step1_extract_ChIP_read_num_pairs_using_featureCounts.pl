#!/usr/bin/perl -w

# Feb 17, 2014
# In this script, inspired by MBD Cell paper,
# I want to use sliding overlaping window approach.

# win size 1kb, moving size 500bp
# From the middle point of DMR (bed region), moving back and forth for 5kb each step is
# 500 bp -9, .., -1, 0, 1, 2,.., 9

#
#ID for featureCounts input is like
# 1:-9
# 1:-8
#1:0
#1:1
#1:9

use strict;
use File::Spec;

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;


print STDERR "no head or head must have first col:  coor or chr \n\n";

my $debug = 0;
my $print_debug = 0;

if ( !$debug ) {
	$print_debug = 0;
}
if($debug){
	print STDERR "debug = 1\n\n";
}

#my $half_interval_length = 2000;
# order for bins: -25,-24,...,-1,0, 1,...25;
my $half_flanking_bin_num = 9;#25; #100; #40;
my $bin_size = 1000;#200;#50;
my $sliding_size = 500; 

my %types = ("bed"=>1, "coor" => 1);

#my $usage = "$0 \n <input_list_DMR> <type: bed or coor> <bam_file>  <bam_label>  <outdir> <output_pre> <CPU_thread>\n\n";
#die $usage unless(@ARGV == 6);
#time /Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/FPKM_based_on_featureCounts/step1_extract_ChIP_read_num_pairs_using_featureCounts.pl coor_nohead.txt coor  bam/ . debug_coor_nohead 1

my $usage = "$0 \n <input_list_DMR> <type: bed or coor> <bam_dir>  <outdir> <output_pre> <CPU_thread>\n\n";
die $usage unless(@ARGV == 6);

my $input	= shift or die;
my $type	= shift or die;
die "type: bed or coor" unless (defined $types{$type});
my $bam_dir	= shift or die;
my $outdir	= shift or die;
my $outpre	= shift or die;

my $CPU_thread = shift or die;

#my $output	= shift or die;
#my $output = File::Spec->catfile($outdir, $outpre . "_half" . $half_interval_length . "bp" . "_bin" . $bin_size ."bp.txt");
my $output = File::Spec->catfile($outdir, $outpre .  "_halfBinN" . $half_flanking_bin_num  . "_bin" . $bin_size ."bp" . "_sliding" . $sliding_size ."bp_PairsN_avg.txt");

die unless (-e $input);
die if( -e $output);
die unless (-d $bam_dir);


my $SAF_file = "";
my $featureCounts_outfile = "";
my $feature_log_file = "";

my ($volume, $directories, $infile_name) =       File::Spec->splitpath( $input );

if ( $infile_name =~ /(\S+)\./) {
	my $pre = $1;
	$SAF_file = File::Spec->catfile($outdir, $pre . "_SAF.txt");
	die if (-e $SAF_file);
	$featureCounts_outfile = File::Spec->catfile($outdir, $pre.  "_featureCounts_pO.txt");
	die if(-e $featureCounts_outfile);
	$feature_log_file  = File::Spec->catfile($outdir, $pre.  "_featureCounts_pO_log.txt");
	die if(-e $feature_log_file);
}

 
Kai_Module::read_list_and_get_sliding_interval ( $input, $type, $SAF_file, $bin_size, $sliding_size, $half_flanking_bin_num );

#featureCounts -T 8 -p -O -F SAF -a TAIR10_WinSize1000bp_sliding500bp_SAF.txt   -o  TAIR10_WinSize1000bp_sliding500bp_featureCounts_for_MBD7_ChIP_pO.txt  [bam]
die unless ( -e $SAF_file);
Kai_Module::run_featureCounts($SAF_file, $featureCounts_outfile, $feature_log_file, $CPU_thread, $bam_dir, $debug);

Kai_Module::extract_featureCounts_results($featureCounts_outfile, $output, $half_flanking_bin_num, $debug);

#Gene_all_coordinate_SAF.txt
#"featureCounts" "-T" "3" "-p" "-O" "-F" "SAF" "-a" "Gene_all_coordinate_SAF.txt" "-o" "all_Gen
#Geneid  Chr     Start   End     Strand  Length  MBD7_5_My..bam

exit;


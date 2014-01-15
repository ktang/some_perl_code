#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $genome_ref = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BowtieIndex/genome";


my $usage = "$0 \n\n <input_fa> <out_pre> <outdir>\n\n";
die $usage unless(@ARGV == 3);

my $input_fa = shift or die "input_fa";
my $pre = shift or die "pre";
my $outdir = shift or die "outdir";

die "input_fa" unless (-e $input_fa);

#####################
# 1 bowtie
#####################

my $p = 8;
my $v = 0;
my $k = 100;

my $sam_file = File::Spec->catfile($outdir, $pre . "_bowtie_v" . $v . "k" . $k . ".sam");
my $log_file = File::Spec->catfile($outdir, $pre . "_log.txt");

if(!$debug){
	die if(-e $sam_file);
}

my $cmd_bowtie = "bowtie -t -f -p $p -v $v -k $k  $genome_ref $input_fa -S $sam_file  >> $log_file  2>&1";

print STDERR $cmd_bowtie, "\n\n";
if(!$debug){
	`$cmd_bowtie`;
}

#####################
# 2 sam2bam
#####################
my $bam_file = File::Spec->catfile($outdir, $pre . "_bowtie_v" . $v . "k" . $k . ".bam");

if(!$debug){
	die unless(-e $sam_file);
	die if(-e $bam_file);
}

my $cmd_bam = "samtools view -b  -S $sam_file > $bam_file";

print STDERR $cmd_bam, "\n\n";
if(!$debug){
	`$cmd_bam`;
}

#################
# sort bam
#################
my $sorted_bam_pre =File::Spec->catfile($outdir, $pre . "_bowtie_v" . $v . "k" . $k . "_sorted");
#File::Spec->catfile($outdir, $pre . "_bowtie2_sorted");
my $sorted_bam = File::Spec->catfile($outdir, $pre . "_bowtie_v" . $v . "k" . $k . "_sorted.bam");
#File::Spec->catfile($outdir, $pre . "_bowtie2_sorted.bam");

if(!$debug){
	die unless(-e $bam_file);
	die if(-e $sorted_bam);
}

my $cmd_sort = "samtools sort $bam_file $sorted_bam_pre ";

print STDERR $cmd_sort, "\n\n";

if(!$debug){
	`$cmd_sort`;
}

#################
# index bam
#################
my $cmd_index = "samtools index $sorted_bam";

print STDERR $cmd_index, "\n\n";
if(!$debug){
	`$cmd_index`;
}

exit;
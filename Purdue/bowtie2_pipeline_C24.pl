#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $bowtie2 = "/Users/tang58/Software/Bowtie/Bowtie2/bowtie2-2.0.0-beta6/bowtie2";
die unless (-e $bowtie2);

my $index_db = "/Users/tang58/DataBase/TAIR_Col0_genome/index/Bowtie2/C24/C24_7chr";

my $ref_fas = 	"/Users/tang58/DataBase/C24/C24_7chr/C24_7chr.fa";
die unless (-e $ref_fas);

my $ref_fai = 	"/Users/tang58/DataBase/C24/C24_7chr/C24_7chr.fa.fai";
die unless(-e $ref_fai);

my $usage = "$0 <indir> <inpre> <outdir>";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
my $pre = shift or die "pre";
my $outdir = shift or die "outdir";

my $fq_in1 = File::Spec->catfile($indir, $pre . "_1.fq");
my $fq_in2 = File::Spec->catfile($indir, $pre . "_2.fq");

die unless (-e $fq_in1);
die unless (-e $fq_in2);

#####################
# 1 bowtie2
#####################

my $sam_file = File::Spec->catfile($outdir, $pre . "_bowtie2.sam");
my $log_file = File::Spec->catfile($outdir, $pre . "_log.txt");


if(!$debug){
	die if(-e $sam_file);
}

my $cmd_bowtie2 = "$bowtie2 -t  -N 1 --ignore-quals -M 3 -I 90 -X 1000 --sensitive --phred64  -x $index_db -1 $fq_in1 -2 $fq_in2 -S $sam_file  >> $log_file  2>&1";

print STDERR $cmd_bowtie2, "\n\n";

if(!$debug){
	`$cmd_bowtie2`;
}

#####################
# 2 sam2bam
#####################
my $bam_file = File::Spec->catfile($outdir, $pre . "_bowtie2.bam");

if(!$debug){
	die unless(-e $sam_file);
	die if(-e $bam_file);
}

my $cmd_bam = "time samtools view -bt $ref_fai $sam_file > $bam_file";

print STDERR $cmd_bam, "\n\n";
if(!$debug){
	`$cmd_bam`;
}

#################
# sort bam
#################
my $sorted_bam_pre = File::Spec->catfile($outdir, $pre . "_bowtie2_sorted");
my $sorted_bam = File::Spec->catfile($outdir, $pre . "_bowtie2_sorted.bam");

if(!$debug){
	die unless(-e $bam_file);
	die if(-e $sorted_bam);
}

my $cmd_sort = "time samtools sort $bam_file $sorted_bam_pre ";

print STDERR $cmd_sort, "\n\n";

if(!$debug){
	`$cmd_sort`;
}

#################
# index bam
#################
my $cmd_index = "time samtools index $sorted_bam";

print STDERR $cmd_index, "\n\n";
if(!$debug){
	`$cmd_index`;
}

#################
# pileup
#################
 my $pileup_file = File::Spec->catfile($outdir, $pre . "_bowtie2.pileup");
 
 if(!$debug){
	die unless(-e $sorted_bam);
	die if(-e $pileup_file);
}
 
 my $cmd_pileup = "samtools pileup -f $ref_fas $sorted_bam | pileup_filter.pl > $pileup_file";
 
print STDERR $cmd_pileup, "\n\n";
if(!$debug){
	`$cmd_pileup`;
}

exit;
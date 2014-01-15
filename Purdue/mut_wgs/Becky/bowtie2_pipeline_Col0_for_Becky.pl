#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

my %Phred_h = ( 33=> "--phred33", 64=>"--phred64");

# my $bowtie2 = "/Users/tang58/Software/Bowtie/Bowtie2/bowtie2-2.0.0-beta6/bowtie2";
my $bowtie2 = "/Users/tang58/Software/Bowtie/Bowtie2/bowtie2-2.0.4/bowtie2";
die "bowtie2" unless (-e $bowtie2);

my $index_db = "/Users/tang58/DataBase/TAIR_Col0_genome/index/Bowtie2/Col0/TAIR10_Col0_7chr";
#die "index_db" unless (-d $index_db);

my $ref_fas = 	"/Users/tang58/DataBase/TAIR_Col0_genome/index/TAIR10_Col0_7chr/TAIR10_Col0_7chr.fa";
die "ref_fas" unless (-e $ref_fas);

my $ref_fai = 	"/Users/tang58/DataBase/TAIR_Col0_genome/index/TAIR10_Col0_7chr/TAIR10_Col0_7chr.fa.fai";
die "ref_fai" unless(-e $ref_fai);

my $usage = "$0 <indir> <inpre> <outdir> <Phred(33/64)>";
die $usage unless(@ARGV == 4);

my $indir = shift or die "indir";
my $pre = shift or die "pre";
my $outdir = shift or die "outdir";

my $phred = shift or die "Phred";

die "Phred = 33 or 64"  unless (defined $Phred_h{ $phred });

my $phred_label = $Phred_h{ $phred };


my $gz1 = File::Spec->catfile($indir, $pre . "_1.fq.gz");
my $gz2 = File::Spec->catfile($indir, $pre . "_2.fq.gz");

if(!$debug){
	die unless (-e $gz1);
	die unless (-e $gz2);
}

#####################
# 1 bowtie2
#####################

my $sam_file = File::Spec->catfile($outdir, $pre . "_bowtie2.sam");
my $log_file = File::Spec->catfile($outdir, $pre . "_log.txt");


if(!$debug){
	die if(-e $sam_file);
}

# --phred64 
my $cmd_bowtie2 = "time $bowtie2 -t -I 0 -X 1000 --sensitive $phred_label  -x $index_db -1 $gz1 -2 $gz2 -S $sam_file  >> $log_file  2>&1";

open (OUT, ">>$log_file") or die;
print OUT "\n", $cmd_bowtie2, "\n\n";
close OUT;

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

open (OUT, ">>$log_file") or die;
print OUT "\n", $cmd_bam, "\n\n";
close OUT;

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
open (OUT, ">>$log_file") or die;
print OUT "\n",$cmd_sort , "\n\n";
close OUT;

print STDERR $cmd_sort, "\n\n";

if(!$debug){
	`$cmd_sort`;
}

#################
# index bam
#################
my $cmd_index = "time samtools index $sorted_bam";

open (OUT, ">>$log_file") or die;
print OUT "\n", $cmd_index, "\n\n";
close OUT;

print STDERR $cmd_index, "\n\n";
if(!$debug){
	`$cmd_index`;
}

exit;
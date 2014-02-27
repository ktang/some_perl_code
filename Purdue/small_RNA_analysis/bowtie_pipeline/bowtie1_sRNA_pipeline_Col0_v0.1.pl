#!/usr/bin/perl -w

# v0.1
# input should be fa.gz
# 

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $genome_ref = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BowtieIndex/genome";


#my $usage = "$0 \n\n <input_fa> <out_pre> <outdir>\n\n";
#die $usage unless(@ARGV == 3);

my $usage = "$0 \n\n <indir> <input_fa_gz_pre> <outdir> <outpre>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die;
my $inpre = shift or die "input_fa";
my $outdir = shift or die "outdir";
my $outpre = shift or die "pre";

#die "input_fa" unless (-e $input_fa);

opendir(DIR, $indir) or die;

my @fa_files = grep /fa\.gz$/, grep /$inpre/ , readdir DIR;

closedir DIR;

my @full_fa_files = map { File::Spec->catfile($indir, $_) } @fa_files;
foreach my $tmp (@full_fa_files){
	die $tmp unless( -e $tmp);
}

my $all_fa_files = join(" ", @full_fa_files);

#####################
# 1 bowtie
#####################

my $p = 2;
my $v = 0;
my $k = 100;

#my $sam_file = File::Spec->catfile($outdir, $pre . "_bowtie_v" . $v . "k" . $k . ".sam");
my $bam_file = File::Spec->catfile($outdir, $outpre . "_bowtie_v" . $v . "k" . $k . ".bam");
my $log_file = File::Spec->catfile($outdir, $outpre . "_log.txt");

if(!$debug){
	#die if(-e $sam_file);
	die if(-e $bam_file);
}

# | perl -lane ' print if( $F[2] ne "*")' | samtools view -b -S - > ./35S_SUC2_Input_131216_bowtie2_Jan14.bam
#my $cmd_bowtie = "zcat $all_fa_files | bowtie -t -f -p $p -v $v -k $k  $genome_ref - -S  >> $log_file  2>&1 | perl -lane ' print if( \$F[2] ne \"*\")' | samtools view -b -S - > $bam_file";
my $cmd_bowtie = "zcat $all_fa_files | bowtie -t -f -p $p -v $v -k $k  $genome_ref - -S  2>> $log_file   | perl -lane ' print if( \$F[2] ne \"*\")' | samtools view -b -S - > $bam_file";

print STDERR $cmd_bowtie, "\n\n";
if(!$debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_bowtie, "\n\n";
	close OUT;
	`$cmd_bowtie`;
}


#################
# sort bam
#################
my $sorted_bam_pre =File::Spec->catfile($outdir, $outpre . "_bowtie_v" . $v . "k" . $k . "_sorted");
#File::Spec->catfile($outdir, $pre . "_bowtie2_sorted");
my $sorted_bam = File::Spec->catfile($outdir, $outpre . "_bowtie_v" . $v . "k" . $k . "_sorted.bam");
#File::Spec->catfile($outdir, $pre . "_bowtie2_sorted.bam");

if(!$debug){
	die unless(-e $bam_file);
	die if(-e $sorted_bam);
}

my $cmd_sort = "samtools sort $bam_file $sorted_bam_pre ";

print STDERR $cmd_sort, "\n\n";

if(!$debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_sort, "\n\n";
	close OUT;
	
	`$cmd_sort`;
}

#################
# index bam
#################
my $cmd_index = "samtools index $sorted_bam";

print STDERR $cmd_index, "\n\n";
if(!$debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_index, "\n\n";
	close OUT;
	
	`$cmd_index`;
}

exit;
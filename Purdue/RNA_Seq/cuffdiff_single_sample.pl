#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;

my $gtf = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf";
die unless (-e $gtf);

my $fa_file = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa";
die unless (-e $fa_file);


if($debug){
#	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <mut_bam> <mut_label> <wt_bam> <wt_label> <outdir> <CPU_num>\n\n";
die $usage unless(@ARGV == 6);

my $bam_mut = shift or die;
my $lab_mut = shift or die;
my $bam_wt  = shift or die;
my $lab_wt  = shift or die;
my $outdir  = shift or die;
my $CPU_num = shift or die;


die unless (-e $bam_mut);
die unless (-e $bam_wt);

die "outdir exists" if (-d $outdir);

my $label = join(",", ($lab_mut, $lab_wt  ));

#my $cmd = "time cufflinks -p $CPU_num --GTF $gtf -o $outdir $bam";
my $cmd = "time cuffdiff -b $fa_file -u $gtf -p $CPU_num -o $outdir -L $label $bam_mut $bam_wt" ;
print STDERR $cmd, "\n\n";
if(!$debug){	
	`$cmd`;
}else{
	print STDERR "OK\n";
}
exit;

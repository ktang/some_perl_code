#!/usr/bin/perl -w


use strict;
use File::Spec;

my $debug = 0;

my $gtf = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf";
die unless (-e $gtf);


if($debug){
#	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <input_bam> <outdir> <CPU_num>\n\n";
die $usage unless(@ARGV == 3);

my $bam = shift or die;
die unless (-e $bam);
my $outdir = shift or die;
die "outdir exists" if (-d $outdir);

my $CPU_num = shift or die;

my $cmd = "time cufflinks -p $CPU_num --GTF $gtf -o $outdir $bam";
print STDERR $cmd, "\n\n";
if(!$debug){	
	`$cmd`;
}else{
	print STDERR "OK\n";
}
exit;

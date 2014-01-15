#!/usr/bin/perl -w

#time tophat2 -p 6
#-G /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf
#-o Fengqiu_WT_thout
# --b2-very-sensitive
#--phred64-quals
#-r 200

#/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome
#CleanData/WT_L1_1.fq.gz CleanData/WT_L1_2.fq.gz 
use strict;
use File::Spec;

my $debug = 0;

my %Phred_h = ( 33=> "", 64=>"--phred64-quals");

my $p = 8;
my $r =20;# 200;
my $mode_flag = "--b2-very-sensitive";


my $gtf = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf";
die unless (-e $gtf);

my $genome = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome";



my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";
die $usage unless(@ARGV == 5);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);
my $f1 = shift or die "f1";
my $f2 = shift or die "f2";

my $outdir = shift or die;
die "outdir" if(-d $outdir);

my $phred = shift or die "Phred";
die "Phred = 33 or 64"  unless (defined $Phred_h{ $phred });
my $phred_label = $Phred_h{ $phred };

my $fq1 = File::Spec->catfile($indir, $f1);
my $fq2 = File::Spec->catfile($indir, $f2);

die unless (-e $fq1);
die unless (-e $fq2);
#time tophat2 -p 6
#-G /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf
#-o Fengqiu_WT_thout
# --b2-very-sensitive
#--phred64-quals
#-r 200

#/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome
#CleanData/WT_L1_1.fq.gz CleanData/WT_L1_2.fq.gz 

my $cmd = "tophat2 -p $p -r $r -G $gtf -o $outdir $mode_flag $phred_label $genome $fq1 $fq2";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
exit;


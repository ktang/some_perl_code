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

#print STDERR (qq/
#	--b2-very-fast
#	--b2-fast
#	--b2-very-sensitive
#/) ;

my %modes_h = (
		"very-fast" 	=> "--b2-very-fast",	
		"fast"		=> "--b2-fast",	
		"sensitive" 	=> "--b2-sensitive",	
		"very-sensitive"=> "--b2-very-sensitive",	
);

my $p = 4;
my $r = 200;

# my $mode_flag = "--b2-very-sensitive";


my $gtf = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf";
die unless (-e $gtf);

my $genome = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome";

my $usage = "$0 \n<indir>  <outdir> <phred_33_64> <mapping_mode>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);
#my $f1 = shift or die "f1";
#my $f2 = shift or die "f2";

my $outdir = shift or die;
die "outdir" if(-d $outdir);

my $phred = shift or die "Phred";
die "Phred = 33 or 64"  unless (defined $Phred_h{ $phred });
my $phred_label = $Phred_h{ $phred };

my $mapping_mode = shift or die;
die "mapping mode" unless (defined $modes_h{$mapping_mode});
 my $mode_flag = $modes_h{$mapping_mode};
 
opendir (INDIR, $indir ) or die;
my @temp = grep /\.gz$/ , readdir INDIR;
closedir INDIR;

die if (@temp != 1);


my $fq1 = File::Spec->catfile($indir, $temp[0]);
die unless (-e $fq1);
#my $fq2 = File::Spec->catfile($indir, $f2);

my $cmd = "tophat2 -p $p -r $r -G $gtf -o $outdir $mode_flag $phred_label $genome $fq1";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
exit;


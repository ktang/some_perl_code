#!/usr/bin/perl -w

# use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/tophat_pair_v0.1.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

print STDERR "change mode_flag p and r in script \n\n";

my $debug = 0;
my $usage = "$0 \n<indir> <all_outdir> < 33 or 64 (phred_33_64) >\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);

my $outdir = shift or die;
die "outdir" unless(-d $outdir);

my $phred = shift or die "phred";

opendir(DIR, $indir) or die "dir";

#my @files = grep /_1.fq.gz/, readdir DIR;

#my @dirs = grep /mRNA_seq/i ,  readdir DIR;
my @dirs = grep /RNA_seq/i ,  readdir DIR;

closedir DIR;

#print STDERR join("\n", @files), "\n\n";
print STDERR join("\n", @dirs), "\n\n";


# my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";

#foreach my $in1(@files){
foreach my $dir (@dirs){
	if($dir =~ /(\S+)_RNA_seq/){
		
		my $pre = $1;
		#my $pre = $dir;
#		my $in2 = $pre . "_l1_2.fq.gz";
		
		my $sub_outdir = File::Spec->catfile($outdir, $pre . "_thout" );
		
		die if (-d $sub_outdir);
		my $sub_indir = File::Spec->catfile($indir, $dir );
		die unless (-d  $sub_indir );
		
		my $cmd = "perl $script $sub_indir  $sub_outdir $phred";
		#"$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}
	}
}

exit;

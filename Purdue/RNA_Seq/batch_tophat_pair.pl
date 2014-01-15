#!/usr/bin/perl -w

# use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/tophat_pair.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

my $debug = 0;
my $usage = "$0 \n<indir> <all_outdir> <phred_33_64>\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);

my $outdir = shift or die;
die "outdir" unless(-d $outdir);

my $phred = shift or die "phred";

opendir(DIR, $indir) or die "dir";

my @files = grep /_1.fq.gz/, readdir DIR;

closedir DIR;

print STDERR join("\n", @files), "\n\n";


# my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";

foreach my $in1(@files){
	if($in1 =~ /(\S+)_l1_1.fq.gz$/){
		
		my $pre = $1;
		
		my $in2 = $pre . "_l1_2.fq.gz";
		
		my $sub_dir = File::Spec->catfile($outdir, $pre . "_thout" );
		
		die if (-d $sub_dir);
		
		my $cmd = "perl $script $indir $in1 $in2 $sub_dir $phred";
		#"$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}
	}
}

exit;
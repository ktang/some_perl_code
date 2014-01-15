#!/usr/bin/perl -w

# use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/cufflinks_single_sample.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

my $debug = 0;
my $usage = "$0 \n<indir> <all_outdir> <CPU_num>\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);

my $outdir = shift or die;
die "outdir" unless(-d $outdir);

my $CPU_num = shift or die;

opendir(DIR, $indir) or die "dir";

my @files = grep /\.bam$/, readdir DIR;

closedir DIR;

print STDERR join("\n", @files), "\n\n";


# my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";

foreach my $file(@files){
	if($file =~ /(\S+)_thout\.bam$/){
		
		my $pre = $1;
		my $sub_dir = File::Spec->catfile($outdir, $pre . "_clout" );
		die if (-d $sub_dir);
		
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		my $cmd = "perl $script $input $sub_dir $CPU_num";
		#"$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}else{
			print STDERR "OK\n";
		}
	}
}

exit;
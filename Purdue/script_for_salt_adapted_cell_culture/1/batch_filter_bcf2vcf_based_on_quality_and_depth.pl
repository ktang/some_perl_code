#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/src/filter_bcf2vcf_based_on_quality_and_depth.pl";
die unless (-e $script);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);


my $indir = shift or die "indir";
my $outdir = shift or die "outdir";

die unless (-d $indir);
die unless (-d $outdir);

opendir (DIR, $indir) or die;
my @files = grep /\.bcf$/ , readdir DIR;
closedir DIR;

print STDERR "files: \n";
print STDERR join("\n", @files ), "\n\n";

foreach my $file (@files){
	my $pre;
	
	if ($file =~ /(\S+)\.raw\.bcf$/){
		$pre = $1;
	}else{
		die $file;
	}
	
	my $input = File::Spec->catfile($indir, $file);
	my $vcf_file = File::Spec->catfile($outdir, $pre . ".flt_Mar28.vcf" );
	my $cmd = "perl $script $input $vcf_file";

	print STDERR $cmd, "\n\n";
#	print STDERR $vcf_cmd, "\n\n\n";
	
	if(!$debug){
		`$cmd`;
	}
}
exit;

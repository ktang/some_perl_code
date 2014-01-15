#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $script = "/Users/tang58/misc/Zhu_Xiaohong/Jun11_start_from_4056loci/src/categorize_split_file_to_get_mutation_type.pl";
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
my @files = grep /\.txt$/ , readdir DIR;
closedir DIR;

print STDERR "files: \n";
print STDERR join("\n", @files ), "\n\n";

foreach my $file (@files){
	my $pre;
	
	if ($file =~ /(\S+)\.txt$/){
		$pre = $1;
	}else{
		die $file;
	}
	
	my $input = File::Spec->catfile($indir, $file);
	my $output = File::Spec->catfile($outdir, $pre . "_with_mut_type.txt" );
	my $cmd = "perl $script $input  $output";

	print STDERR $cmd, "\n\n";
#	print STDERR $vcf_cmd, "\n\n\n";
	
	if(!$debug){
		`$cmd`;
	}
}
exit;

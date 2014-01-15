#!/usr/bin/perl -w

use strict;
use File::Spec;

my $script = "/Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/mpileup/513loci_downstream/self_percentage_cutoff_Jun13/src/filter_SelfPer_input_col_number_and_cutoff.pl";
die unless (-e $script);

#my $usage = "\n$0 <input> <cutoff> <col_number> STDOUT\n\n";

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> <cutoff> <col_number>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die "indir";
my $outdir = shift or die "outdir";
my $cutoff = shift or die;
my $col_num = shift or die;

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
	my $output = File::Spec->catfile($outdir, $pre . "_SelfPer" . $cutoff . ".txt" );
	die if(-e $output);
	my $cmd = "perl $script $input $cutoff $col_num > $output";

	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}
exit;
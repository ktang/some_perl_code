#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0 ;

my $script = "/Users/tang58/scripts_all/perl_code/count_overlap_for_Venn/count_overlap_info_for_Venn.pl";
die unless (-e $script);

print STDERR "\n\n script used is\n";
print STDERR $script, "\n\n";

my $usage = "$0 \n<indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;
die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir);
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

#<file_name> <indir> <outdir> <output_name>
foreach my $file( @files ){
	my $output;
	my $input = File::Spec->catfile($indir, $file);
	die unless (-e $input);
	
	if($file =~ /(\S+)\.txt$/){
		$output = File::Spec->catfile($outdir, $1 . "_num_Venn.txt");
		die if (-e $output);
	}else{
		die $file;
	}
	
	my $cmd = "perl $script $input $output";
#	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

exit;

#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/coor2bed/more_bed2bed_5col.pl";

die unless (-e $script);

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <indir> <outdir>\n\n";
#die $usage unless(@ARGV == 2);

my $usage = "$0 \n <indir> <outdir> <chr N>\n\n";
die $usage unless(@ARGV == 3);


my $indir = shift or die;
my $outdir = shift or die;

my $chr_N = shift or die;
die $usage unless ($chr_N eq "chr" or $chr_N eq "N" );


die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir) or die "cannot open $indir: $!";
my @files = grep /txt$/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	if ( $file =~ /(\S+)\.txt$/) {
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $output = File::Spec->catfile($outdir, $pre . "_coordinate_bed.txt");
		my $cmd = "perl $script $input $output $chr_N";
		print STDERR $cmd, "\n\n";
		if (!$debug) {
			`$cmd`;
		}
		
	}
	
}

exit;

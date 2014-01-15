#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 0;
my $min_quality = 700;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/BAM_SAM_operation/filter_bam_with_insert_size_for_distant_pair_v0.1.pl";
die unless (-e $script );

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> \n\n";
die $usage unless(@ARGV == 2);

#my $input = shift or die;
my $indir  = shift or die;
my $outdir = shift or die;

die unless (-d $indir) ;
die unless (-d $outdir);

opendir (INDIR, $indir) or die;
my @files = grep /\.bam$/, readdir INDIR;
closedir INDIR;

foreach my $file (@files){
#	if( $file =~ /(\S+)\.bam$/){
	if( $file =~ /(\S+)_Insert700/){
		my $output = File::Spec->catfile($outdir , $1 . "_insertSize" . $min_quality . "_only.bam");
		die if(-e $output);
		my $input = File::Spec->catfile ( $indir, $file );
		die unless (-e $input);
		my $cmd = "time perl $script $input $output $min_quality";
		print STDERR $cmd, "\n\n";
		if(!$debug){
			`$cmd`;
		}
	}else{
		die $file;
	}
}

exit;

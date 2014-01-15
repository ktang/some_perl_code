#!/usr/bin/perl -w


my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/single_read/tophat_single_end.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

my %modes_h = (
		"very-fast" 	=> "--b2-very-fast",	
		"fast"		=> "--b2-fast",	
		"sensitive" 	=> "--b2-sensitive",	
		"very-sensitive"=> "--b2-very-sensitive",	
);

my $debug = 0;
my $usage = "$0 \n<indir> <all_outdir> <phred_33_64> <mapping_mode>\n\n";
die $usage unless(@ARGV == 4);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);

my $outdir = shift or die;
die "outdir" unless(-d $outdir);

my $phred = shift or die "phred";

my $mapping_mode = shift or die;

die (qq/
	--b2-very-fast
	--b2-fast
	--b2-very-sensitive
/) unless (defined $modes_h{$mapping_mode});

opendir(DIR, $indir) or die "dir";
my @dirs = grep /RNA_seq/ ,  readdir DIR;
closedir DIR;

#print STDERR join("\n", @files), "\n\n";
print STDERR join("\n", @dirs), "\n\n";
foreach my $dir (@dirs){
	if($dir =~ /(\S+)_RNA_seq/){
		
		my $pre = $1;
		
#		my $in2 = $pre . "_l1_2.fq.gz";
		
		my $sub_outdir = File::Spec->catfile($outdir, $pre . "_thout" );
		
		die if (-d $sub_outdir);
		my $sub_indir = File::Spec->catfile($indir, $dir );
		die unless (-d  $sub_indir );
		
		my $cmd = "perl $script $sub_indir  $sub_outdir $phred $mapping_mode";
		#"$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}
	}
}

exit;
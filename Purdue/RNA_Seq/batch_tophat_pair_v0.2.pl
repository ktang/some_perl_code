#!/usr/bin/perl -w

#v0.2
# input mapping mode

# Preset options in --end-to-end mode (local alignment is not used in TopHat2)
#    --b2-very-fast
#     --b2-fast
#      --b2-sensitive
#     --b2-very-sensitive


# use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

print STDERR "incomplete\n\n";
exit;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/tophat_pair_v0.2.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

print STDERR (qq/
	--b2-very-fast
	--b2-fast
	--b2-sensitive
	--b2-very-sensitive
/) ;

my %modes_h = (
		"very-fast" 	=> "--b2-very-fast",	
		"fast"		=> "--b2-fast",	
		"sensitive" 	=> "--b2-sensitive",	
		"very-sensitive"=> "--b2-very-sensitive",	
);

my $debug = 0;
my $usage = "$0 \n<indir> <all_outdir> < 33 or 64 (phred_33_64) > <mode>\n\n";
die $usage unless(@ARGV == 3);

my $indir = shift or die "indir";
die "indir" unless (-d $indir);

my $outdir = shift or die;
die "outdir" unless(-d $outdir);

my $phred = shift or die "phred";
die "phred = 33 or 64" unless ($phred == 33 or $phred == 64);

my $mode = shift or die;
die (qq/
	--b2-very-fast
	--b2-fast
	--b2-sensitive
	--b2-very-sensitive
/) unless (defined $modes_h{$mode});

opendir(DIR, $indir) or die "dir";

#my @files = grep /_1.fq.gz/, readdir DIR;

my @dirs = grep /RNA_seq/ ,  readdir DIR;

closedir DIR;

#print STDERR join("\n", @files), "\n\n";
print STDERR join("\n", @dirs), "\n\n";


# my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";

#foreach my $in1(@files){
foreach my $dir (@dirs){
	if($dir =~ /(\S+)_RNA_seq/){
		
		my $pre = $1;
		
#		my $in2 = $pre . "_l1_2.fq.gz";
		
		my $sub_outdir = File::Spec->catfile($outdir, $pre . "_thout" );
		
		die if (-d $sub_outdir);
		my $sub_indir = File::Spec->catfile($indir, $dir );
		die unless (-d  $sub_indir );
		
		my $cmd = "perl $script $sub_indir  $sub_outdir $phred $mode";
		#"$dump -Q 64 --defline-seq \"@\\\$sn\" --defline-qual \"+\" -E $input";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}
	}
}

exit;
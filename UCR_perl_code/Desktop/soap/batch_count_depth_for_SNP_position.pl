#!/usr/bin/perl -w
# use single_count_depth_for_SNP_position.pl to batch work

use strict;

my $usage = "$0 <indir> <out_dir>";

die $usage unless(@ARGV == 2);

my ($indir, $outdir) = @ARGV[0..1];

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my @files = grep {/\.soapout$/} readdir INDIR;

foreach my $file(@files){
     if($file =~ /(\S+)\.soapout$/){
		 my $pre = $1;
		 my $output = $pre . "_SNP_depth";
		 my $cmd = "single_count_depth_for_SNP_position.pl $indir/$file $outdir/$output";
		 
		 print $cmd, "\n";
		 `$cmd`;
	 }
}

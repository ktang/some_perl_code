#!/usr/bin/perl -w
# batch generate gff file from preprocessed file
use strict;

my $usage = "$0 <id_pos_file> <preprocessed_dir> <gff_dir>";
die $usage unless(@ARGV == 3);

my $indir = $ARGV[1];
my $outdir = $ARGV[2];

opendir (INDIR, $ARGV[1]) or die "Cannot open dir $indir:$!";
my @files = grep {/d\.txt$/} readdir INDIR;
foreach my $file(@files){
     if($file =~ /(\S+)\.txt$/){
		 my $pre = $1;
		 my $output = $pre . ".gff";
		 my $cmd = "perl ~/Desktop/perl_code/NimbleGen/generate_gff.pl $ARGV[0] $indir/$file $outdir/$output";
		 
		 print $cmd, "\n";
		 my $msg = `$cmd`;
		 print $msg, "\n";
	 }
}

print "\a";
exit;
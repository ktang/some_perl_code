#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
use File::Spec;
use strict;

my $script = "/Volumes/Macintosh_HD_2/Histone_modification/GSE28398/SRA/map_ChIP_seq_SOLiD_in_dir.pl";
die unless (-e $script);

my $debug = 1;
my $usage = "$0 <indir>";
die $usage unless(@ARGV == 1);

my $indir = shift or die;
die unless (-d $indir);

opendir ( DIR, $indir ) or die;

my @dirs = grep /^GSM70/, readdir DIR;

foreach my $dir ( @dirs ){
	my $label;
	
	if($dir =~ /^GSM70\d+_(\S+)$/){
		$label = $1;
	}else{
		die $dir;
	}
	
	
	my $cmd = "time perl $script $dir $label";
	print STDERR $cmd, "\n\n";

	if(!$debug){
		`$cmd`;
	}
}

exit;

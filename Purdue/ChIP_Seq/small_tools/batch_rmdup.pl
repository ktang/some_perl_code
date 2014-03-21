#!/usr/bin/perl -w

# this script take a dir as input, in the dir, there only
# sorted paired-end bam files,
# just rmdup

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "\n$0 \n\n<indir> \n\n";
die $usage unless(@ARGV == 1);

my $indir     = shift or die;

opendir(DIR, $indir) or die;
#my @files = grep /sort/, grep /\.bam$/, readdir DIR;
my @files =  grep /\.bam$/, readdir DIR;
closedir DIR;

foreach my $file (@files) {
	if( $file =~ /(\S+)\.bam$/){
		my $outpre = $1;
		my $log_file = File::Spec->catfile($indir, $outpre . "_rmdup_log.txt");
		die "log file exist\n\n" if(-e $log_file);
		my $sorted_bam  = File::Spec->catfile($indir, $file);
		die unless (-e $sorted_bam);
		my $bam_rmdup = File::Spec->catfile($indir, $outpre . "_rmdup.bam" );

		#######
		# rmdup
		########

		if(!$debug){
			die unless(-e $sorted_bam);
			die if(-e $bam_rmdup);
		}
		my $cmd_rmdup = "samtools rmdup $sorted_bam $bam_rmdup 2>> $log_file";
		print STDERR $cmd_rmdup, "\n\n";
		unless($debug){
			open(OUT, ">>$log_file") or die;
			print OUT $cmd_rmdup, "\n";
			close OUT;
			`$cmd_rmdup`;
		}

		######
		# index
		##########

		my $cmd_index = "samtools index $bam_rmdup";
		print STDERR $cmd_index, "\n\n";
		unless($debug){
			open(OUT, ">>$log_file") or die;
			print OUT "\n\n", $cmd_index, "\n\n";
			close OUT;
			`$cmd_index`;
		}
		
	}else{
		die $file;
	}
}
exit;
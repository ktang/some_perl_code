#!/usr/bin/perl -w
 

use strict;

my $usage = "$0 <dir>";
die $usage unless(@ARGV == 1 );

my $indir = $ARGV[0];
#my ($indir, $db) = @ARGV[0..1];
my $db = "/Users/tang58/DataBase/TAIR_Col0_genome/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.index";
opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my @files = grep {/\.fastq$/} readdir INDIR;



foreach my $file(@files){
	if($file =~/(\S+).fastq$/ ){
		my $pre = $1;
		my $unmapped = $pre."_unmapped_to_genome.fasta";
		if (-e "$indir/$unmapped"){die "$unmapped exists!!!"}
		my $cmd = "soap -a $indir/$file -D $db -n 1 -r 2 -v 2 -M 4 -o /dev/null -u $indir/$unmapped 2> /dev/null";
		print $cmd, "\n\n";
		system("$cmd");
	}
}


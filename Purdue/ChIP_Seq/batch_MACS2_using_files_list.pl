#!/usr/bin/perl -w

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $genome_size = "1.2e8";
my $keep_dup = "all";

#my $usage = "\n$0 \n\n <treatment_ChIP> <control_Input> <outdir> <outpre> \n\n";
#die $usage unless(@ARGV == 4);

my $usage = "\n$0 \n\n <input_list> <indir> <outdir> \n\n";
die $usage unless(@ARGV == 3);


my $inlist = shift or die;
my $indir  = shift or die;

#my $tfile     = shift or die;
#my $cfile     = shift or die;
my $outdir    = shift or die;
#my $outpre    = shift or die;

die unless (-d $indir );
die unless (-d $outdir );

open(IN, $inlist) or die;
while (<IN>) {
	chomp;
	my @a = split "\t";
	next if ($a[0] =~ /Treatment/ );
	
	my ($tfile , $cfile) = map {File::Spec->catfile ($indir, $_)} @a[0..1];
	
	my $outpre = $a[2];
	
	my $log_file = File::Spec->catfile($outdir, $outpre . "_MACS2_log.txt");

	die unless (-e $tfile );
	die unless (-e $cfile );
	die if (-e $log_file );
	
	die "log file exist\n\n" if(-e $log_file);

	my $cmd = "macs2 callpeak -t $tfile -c $cfile -g $genome_size --keep-dup $keep_dup --outdir $outdir -n $outpre 2>> $log_file";

	print STDERR $cmd, "\n\n";
	unless($debug){
		open(OUT, ">>$log_file") or die;
		print OUT $cmd, "\n";
		close OUT;
		`$cmd`;
	}
}

close IN;

exit;

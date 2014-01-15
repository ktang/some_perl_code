#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

my $snap_shot_method = "snapshotmainViewWithLabels"; 
#snapshot
#snapshotwholeFrame
#snapshotmainView
#snapshotmainViewWithLabels
#snapshotslicedViewWithLabels

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_list> <output>  <snapshot_outdir>  <snapshot_pre> <up_bp> <down_bp> <sleep_second>\n\n";
die $usage unless(@ARGV == 7);

my $input = shift or die;
my $output = shift or die;

my $snapshot_outdir = shift or die "snapshot_outdir";
my $snapshot_pre    = shift or die "snapshot_pre";

#(my $up_bp          = shift) or ( die "up_bp" if ($up_bp != 0));
#(my $down_bp		 = shift) or ( die "down_bp" if ($down_bp != 0));

my $up_bp          = shift;
die "up_bp" unless ( defined $up_bp );
my $down_bp		 = shift;
die "down_bp" unless  (defined $down_bp );

my $sleep_second = shift or die "sleep_second";
$sleep_second *= 1000;

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $h = <IN>;

my $i = 0;

while (<IN>){
	chomp;
	$i++;
	my @a = split "\t";
	my ($chr, $start, $end ) = @a[0..2];
	my $label = join("_",($chr, $start, $end )  );
	my $file_name = File::Spec->catfile($snapshot_outdir, $snapshot_pre . "_$i" . "_$label" . ".png");
	my $start_up = ( ($start - 1 - $up_bp) >=0 ) ? ($start - 1 - $up_bp)  : 0;
	my $region = $chr . ":" . $start_up . "-" . ( $end - 1 + $down_bp);
	print OUT "goto $region\n";
	print OUT "refresh\n";
	print OUT "sleep $sleep_second\n";	
	print OUT "$snap_shot_method $file_name\n"; 
}



close(IN);
close(OUT);

exit;
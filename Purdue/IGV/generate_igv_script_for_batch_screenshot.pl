#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

#v0.1
# use mid-point then cut exactly XXXbp around the mid-point

use strict;
use File::Spec;

my $debug = 0;

#snapshotDirectory /Users/tang58/try/Learn/IGV

#goto chr1:8666-8710
#collapse
#snapshot

#goto chr1:270494-270622
#collapse
#snapshot


#my $snap_shot_method = "snapshotmainViewWithLabels"; 
#snapshot
#snapshotwholeFrame
#snapshotmainView
#snapshotmainViewWithLabels
#snapshotslicedViewWithLabels

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <input_list> <snapshot_outdir>  <snapshot_pre> <exact_bp> <sleep_second>  <output> \n\n";
#die $usage unless(@ARGV == 6);

my $usage = "$0 \n <input_list> <snapshot_outdir>  <exact_bp> <output_igb_script>\n\n";
die $usage unless(@ARGV == 4);


my $input = shift or die;

my $snapshot_outdir = shift or die "snapshot_outdir";
#my $snapshot_pre    = shift or die "snapshot_pre";

#(my $up_bp          = shift) or ( die "up_bp" if ($up_bp != 0));
#(my $down_bp		 = shift) or ( die "down_bp" if ($down_bp != 0));

#my $up_bp          = shift;
#die "up_bp" unless ( defined $up_bp );
#my $down_bp		 = shift;
#die "down_bp" unless  (defined $down_bp );

my $bp = shift or die;

#my $sleep_second = shift or die "sleep_second";
#$sleep_second *= 1000;

my $output = shift or die;


die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

print OUT "snapshotDirectory $snapshot_outdir\n";

my $h = <IN>;

my $i = 0;

my $half = int ( $bp / 2 );

while (<IN>){
	chomp;
	$i++;
	my @a = split "\t";
	my ($chr, $start, $end ) = @a[0..2];
	my $label = join("_",($chr, $start, $end )  );
#	my $file_name = File::Spec->catfile($snapshot_outdir, $snapshot_pre . "_$i" . "_$label" . ".png");
#	my $start_up = ( ($start - 1 - $up_bp) >=0 ) ? ($start - 1 - $up_bp)  : 0;
#	my $region = $chr . ":" . $start_up . "-" . ( $end - 1 + $down_bp);
	
	my $mid_point = int( ($start + $end) / 2 );
	
	my $cut_start = 0;
	if( $mid_point - $half - 1 >= 0) {
		$cut_start = 	$mid_point - $half - 1;
	}
	
	my $cut_end = $mid_point + $half - 1;
	
	my $region = $chr . ":" . $cut_start . "-" . $cut_end;
	
	print OUT "goto $region\n";
#	print OUT "refresh\n";
#	print OUT "sleep $sleep_second\n";
	print OUT "collapse\n";
	print OUT "snapshot\n"; 
}



close(IN);
close(OUT);

exit;
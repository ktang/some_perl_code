#!/usr/bin/perl -w

#v0.1
# use mid-point then cut exactly XXXbp around the mid-point

#v0.2 input is fasta head like
# >AT1G08060.1 | Symbols: MOM, MOM1 | ATP-dependent helicase family protein | chr1:2501742-2511084 REVERSE LENGTH=9343

use strict;
use File::Spec;

my $debug = 0;

my $up = 2000;
my $down = 500;

my $snap_shot_method = "snapshotmainViewWithLabels"; 
#snapshot
#snapshotwholeFrame
#snapshotmainView
#snapshotmainViewWithLabels
#snapshotslicedViewWithLabels

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_file> <snapshot_outdir>  <sleep_second>  <output> \n\n";
die $usage unless(@ARGV == 4);

my $input = shift or die;

my $snapshot_outdir = shift or die "snapshot_outdir";
#my $snapshot_pre    = shift or die "snapshot_pre";

my $sleep_second = shift or die "sleep_second";
$sleep_second *= 1000;

my $output = shift or die;


die unless (-e $input);
die if( -e $output);



open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $i = 0;
while (<IN>){
	chomp;
	$i++;
	
	my @a  = split /\|/;
	
	my $coordinate ;
	my $strand;
	my $symbol = "unknown";
	my $id;
	
	if(/(AT\dG\d\d\d\d\d)\.\d/){
		$id = $1;
	}else{
		die "id $_";
	}
	
	if (/(chr\d:\d+-\d+)\s+(\S+)\s+LENGTH/){
		$coordinate = $1;
		$strand = $2;
	}else{
		die "region, $_";
	}
	
	die unless ( $strand eq "FORWARD" or $strand eq "REVERSE" );
	
	if ( $a[1] =~ /Symbols:\s+(\w+),*\s+/){
		$symbol = $1;
	}else{
		print STDERR $id , "\n";
	}
	
	
	
	my ($chr, $first, $second) ;
	
	if ($coordinate =~ /(chr\d):(\d+)-(\d+)/){
		($chr, $first, $second) = ($1, $2, $3);
	}else{
		die $coordinate;
	}
	
	my $region;
	if( $strand eq "FORWARD" ){
		$region = $chr . ":" . ($first - $up)   . "-" . ($second + $down)  ;
	}else{
		$region = $chr . ":" . ($first - $down) . "-" . ($second + $up)  ;
	}
	
#	my $file_name = File::Spec->catfile($snapshot_outdir, $snapshot_pre . "_$i" . "" . ".png");
	my $file_name = File::Spec->catfile($snapshot_outdir, $i . "_" . $id . "_" .  $symbol . ".png");

	print OUT "goto $region\n";
#	print OUT "refresh\n";
	print OUT "sleep $sleep_second\n";	
	print OUT "$snap_shot_method $file_name\n";
	
	if($debug){
		if ($i >= 5){exit}
	#	print STDERR join("\n", ($i, $coordinate, $strand, $symbol, $chr, $first, $second) ) , "\n\n";
	#	print STDERR join("\n", $i, @a), "\n\n\n";
		
	}
}



close(IN);
close(OUT);

exit;
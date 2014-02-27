#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
my $print = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0\n <head: y or n> <depth_cutoff> <gap_cutoff> \n\n";
die $usage unless(@ARGV == 3);

my $head_flag    = shift or die;
my $depth_cutoff = shift or die;
my $gap_cutoff   = shift;
die "gap_cutoff" unless ( $gap_cutoff >= 0);


die $usage unless ( $head_flag eq "y" or  $head_flag eq "n"  );

my ($last_chr, $last_start, $last_end) = ("chr0", -1 , -1);

my $last_pos = -1;

if($head_flag eq "y" ){
	my $h = <>;
	print $h;
}

#while(<IN>){
while(<>){
	chomp;
	my ($chr, $start_l, $end_l, $dep) = split "\t";
	
	next unless ($dep >= $depth_cutoff);
	
	if($chr ne $last_chr ){
		if( $last_chr ne "chr0"){
			print join("\t",($last_chr, $last_start, $last_end)), "\n";
		}
		
		#($last_chr, $last_start, $last_end)  =  ($chr, $pos, $pos);
		($last_chr, $last_start, $last_end)  =  ($chr,  $start_l, $end_l ) if ($dep >= $depth_cutoff );
		  
	}else{
		#if($pos == $last_pos + 1){
		#	$last_end = $pos;
		#}else{
		#	print join("\t",($last_chr, $last_start, $last_end)), "\n";
		#	 ($last_chr, $last_start, $last_end)  =  ($chr, $pos, $pos);
		#}
		
		if ($start_l - $last_end > $gap_cutoff) {
			print join("\t",($last_chr, $last_start, $last_end)), "\n";
			($last_chr, $last_start, $last_end)  =  ($chr,  $start_l, $end_l )
		}else{
			$last_end = $end_l;
		}
		
		
	}
	#$last_pos  = $pos;
}

print join("\t",($last_chr, $last_start, $last_end)), "\n";

#close(IN);

exit;

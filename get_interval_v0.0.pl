#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;
my $print = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#my $usage = "$0 <input>";
#die $usage unless(@ARGV == 1);
#my $input = shift or die;
#die unless (-e $input);
#open(IN, $input) or die;

my $usage = "$0\n <head: y or n> \n\n";
die $usage unless(@ARGV == 1);

my $head_flag = shift or die;

die $usage unless ( $head_flag eq "y" or  $head_flag eq "n"  );

my ($last_chr, $last_start, $last_end) = ("chr0", -1 , -1);

my $last_pos = -1;

if($head_flag eq "y" ){
	my $h = <>;
}

#while(<IN>){
while(<>){
	chomp;
	my ($chr, $pos) = split "\t";
	
	if($chr ne $last_chr ){
		if( $last_chr ne "chr0"){
			print join("\t",($last_chr, $last_start, $last_end)), "\n";
		}
		
		($last_chr, $last_start, $last_end)  =  ($chr, $pos, $pos);
		  
	}else{
		if($pos == $last_pos + 1){
			$last_end = $pos;
		}else{
			print join("\t",($last_chr, $last_start, $last_end)), "\n";
			 ($last_chr, $last_start, $last_end)  =  ($chr, $pos, $pos);
		}
		
		
	}
	
	
	$last_pos  = $pos;
}

#close(IN);

exit;
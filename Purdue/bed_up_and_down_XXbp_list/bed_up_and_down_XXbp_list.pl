#!/usr/bin/perl -w

#bed_up_and_down_XXbp_list.pl
#ReadMe
# input a bed_like list,
# output two list, one is up XXXbp region of the bed_list
# the other is dwonstream XXXbp region of the bed_list


use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <outdir> <XXbp>\n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $outdir = shift or die;
my $bp_num = shift or die;

die unless (-e $input);

my ($volume,$directories,$file) = File::Spec->splitpath( $input );

my ($output_up, $output_down);
if ($file =~ /(\S+)\.txt$/) {
	$output_up   = File::Spec->catfile($outdir, $1 . "_up"   . $bp_num . "bp_region.txt");
	$output_down = File::Spec->catfile($outdir, $1 . "_down" . $bp_num . "bp_region.txt");
}


open(IN, $input) or die "cannot open $input: $!";

die if(-e $output_up);
die if(-e $output_down);
#open(OUT, ">$output") or die "cannot open $output: $!";
open(UP, ">$output_up") or die "cannot open $output_up: $!";
open(DOWN, ">$output_down") or die "cannot open $output_down: $!";

print UP join("\t", ("chr", "start", "end", "origianl_loci")), "\n";
print DOWN join("\t", ("chr", "start", "end", "origianl_loci")), "\n";

while (<IN>) {
	chomp;
	my @a = split "\t";
	next if($a[0] =~ /^chr$/i);
	my ($chr, $start,$end) = @a[0..2];
	$chr = simple_chr($chr);
	
	print UP join("\t", ($chr, (( $start-$bp_num >=1 )?($start-$bp_num ):(1)), ( ($start-1>=1)?($start-1):(1)), "$chr:$start-$end")), "\n";
	print DOWN join("\t", ($chr, $end+1, $end+$bp_num, "$chr:$start-$end")), "\n";
}


close(IN);
close(UP);
close(DOWN);

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
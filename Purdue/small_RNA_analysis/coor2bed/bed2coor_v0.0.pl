#!/usr/bin/perl -w

#v0.0
#input  is just
#chr	start	end
#chr1	34	56
# output file first col is like 1:34-56
# and omit Pt and Mt

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $h = <IN>;
chomp ($h);
my @h_a = split "\t", $h;
print OUT join("\t", ( "coor", @h_a[3..$#h_a]) ) ,"\n";
#print OUT join("\t", ("chr", "start", "end") ), "\n";

#my ($chr, $s, $e) ;

while (<IN>) {
	chomp;
	my @a=split "\t";
	next if ($a[0] =~ /[PMt]/ or $a[0] eq "chrC");
	my $chr = simple_chr($a[0]);

	print OUT join ("\t", ("$chr:$a[1]-$a[2]", @a[3..$#a]) ),"\n";
	
}


close(IN);
close(OUT);

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
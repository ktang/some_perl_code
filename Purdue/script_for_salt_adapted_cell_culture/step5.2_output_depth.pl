#!/usr/bin/perl -w

# step5.2_output_depth.pl
# As I found that in the Batch1, the depth is quite small so that using dep=8 as cutoff,
# All B1 samples seems strange.
#


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
#use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;


my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}


my $usage = "$0 \n <input> <outpu> <sample_num> \n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $output = shift or die;
#my $filted_file = shift or die;
#my $outpre = shift or die;
my $sample_num  = shift or die;
#my $phred_33_64 = shift or die;
 
#die "Phred_33_64" unless ($phred_33_64 == 33 or $phred_33_64 == 64);
#my $output = $outpre . "_dep" . $min_dep . "_per" . $min_per . ".txt";
die unless (-e $input);
die if( -e $output);
#die if( -e $filted_file);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";
#open( POL, ">$filted_file") or die;
my $head = <IN>;
print OUT $head;
#print POL $head;


while (<IN>) {
	chomp;
	my @a = split "\t";
	
#	print OUT join ("\t", @a[0..2]), "\t";
	#my $ref_base = $a[2];
#	my $flag = 0; # if flag = 0; DO NOT output the line; elsif flag = 1 , output it
#	my $benchmark = undef;
	print OUT join("\t", @a[0..2]) , "\t";
 	for my $n ( 1..$sample_num ){
		my $this = $a[$n + 2];
		#my @paires = split ";", $this;
	#DEP=14;TYPE=DEL;ALT=DEL;PER=7.14
	#DEP=17;TYPE=SNP,DEL,INS;ALT=A,DEL,INS;PER=5.88,5.88,11.76
		if ($this =~ /DEP=(\d+)/) {
			print OUT $1;
		}else{
			die $this;
		}
		if ( $n != $sample_num) {
			print OUT "\t";
		}else{
			print OUT "\n";
		}
	}
}
close(IN);
close(OUT);

exit;
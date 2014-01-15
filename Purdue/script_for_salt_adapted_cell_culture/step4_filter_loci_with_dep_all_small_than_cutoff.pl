#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


#BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
#use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;


my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#my $min_base_qual = 20;

my $min_dep = 8;


my $usage = "$0 \n <input> <output> <sample_num> \n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $output = shift or die;
my $sample_num  = shift or die;
#my $phred_33_64 = shift or die;
 
#die "Phred_33_64" unless ($phred_33_64 == 33 or $phred_33_64 == 64);

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $head = <IN>;
print OUT $head;

while (<IN>) {
	chomp;
	my @a = split "\t";
#	print OUT join("\t", (@a[0..2])), "\t";
	my $print_flag = 0;
	for my $n ( 1..$sample_num ){
		my $this_dep = $a[ 3* $n];
		if ($this_dep >= $min_dep) {
			$print_flag = 1;
			last;
		}
	}
	
	if ( $print_flag ) {
		print OUT $_, "\n";
	}
}
close(IN);
close(OUT);

exit;
#!/usr/bin/perl -w

use strict;

my $usage = "$0 <dir>";
die $usage unless (@ARGV == 1);

my $dir = $ARGV[0];

opendir(DIR, $dir);

my @files = grep /\.soapout$/, readdir DIR;

foreach my $file(@files){
	my $cmd = "time /Users/tang58/scripts_all/perl_code/Purdue/Ros3_binding_project/Liu_suggestion/count_soap_genome.pl $dir $file";
	
	print STDERR "$cmd\n\n\n";
	`$cmd`;
}
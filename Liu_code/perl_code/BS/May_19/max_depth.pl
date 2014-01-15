#!/usr/bin/perl -w
use strict;


print "Max Depth:\n";
foreach my $file(@ARGV){
	open(IN, $file) or die "Can't open $file:$!";
    my $max = 0;
    while(<IN>){
	    next if (/SampleID/);
	    my @temp = split /\t/;
	    if($max < $temp[2]){
		    $max = $temp[2];
	    }
	}
	print $file, ": ", $max, "\n";
}

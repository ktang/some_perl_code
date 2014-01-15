#!/usr/bin/perl -w
# from soap output, get coverage
use strict;

my $usage = "$0 <soap_out> [start, end]";
die $usage unless(@ARGV >= 1);
my $soap = $ARGV[0];
my ($start, $end);
if(@ARGV >= 2){
	$start = $ARGV[1];
}
if(@ARGV >= 3){
	$end = $ARGV[2];
}

open(IN, $soap) or die "Can't open $soap: $!";
my (@cover, @neg_cover);
my ($min, $max) = (100000000000000000, 1);
my ($neg_min, $neg_max) = (10000000000000000, 1);
while(<IN>){
	my ($id, $seq, $qual, $hits, $ab, $len, $strand, $chr, $loc, $type) = split /\t/;
    my $end_pos = $loc+$len-1;
	if($strand eq '+'){
		
		if($loc < $min){
			$min = $loc;
		}
		if($end_pos > $max){
			$max = $end_pos;
		}
		foreach my $i($loc..$end_pos){
			$cover[$i]++;
		}
	}elsif($strand eq '-'){
		if($loc < $neg_min){
			$neg_min = $loc;
		}
		if($end_pos > $neg_max){
			$neg_max = $end_pos;
		}
		foreach my $i($loc..$end_pos){
			$neg_cover[$i]++;
		}
	}else{
		die "strand $strand does not match pattern";
	}
}
if(!defined $start){
	$start = $min;
}
if(!defined $end){
	$end = $max;
}
#print "positive strand\n";
foreach my $i($start..$end){
	if(defined $cover[$i]){
	    print $i, "\t", $cover[$i], "\n";
	}else{
		print $i, "\t", "0", "\n";
	}
}
#print "negative strand\n";

foreach my $i($start..$end){
	if(defined $neg_cover[$i]){
	    print $i, "\t-", $neg_cover[$i], "\n";
	}else{
	#	print $i, "\t", "0", "\n";
	}
}


#!/usr/bin/perl -w
# from soap output, count number of small RNAs that map to luc region
# added 8/24/2011: give information on all matched sequences (both positive and negative strands

use strict;

my $usage = "$0 <soap_out> [soap_out 2 ....]";
die $usage unless(@ARGV >= 1);
my ($start, $end) = (745, 2397);

my %seqInfo;
my %seqCnts;
foreach my $soap(@ARGV){

open(IN, $soap) or die "Can't open $soap: $!";
my %reads;
my ($num_pos, $num_neg) = (0,0);
while(<IN>){
	my ($id, $seq, $qual, $hits, $ab, $len, $strand, $chr, $loc, $type) = split /\t/;
    my $end_pos = $loc+$len-1;
	if($loc > $end || $end_pos < $start){
		next;
	}

	$reads{$len}++;
	if($strand eq "-"){
		$seq =~ tr/ACGTacgt/TGCAtgca/;
	}

	#my $key = $seq;
	if(!defined $seqInfo{$seq}){
		$seqInfo{$seq} = [$len, $strand, $loc, $end_pos];
	}else{
		if($seqInfo{$seq}->[0] != $len || $seqInfo{$seq}->[1] ne $strand ||
			$seqInfo{$seq}->[2] != $loc || $seqInfo{$seq}->[3] != $end_pos){
			die "seq info $seq does not match existing info", "current info:",
			"len=$len, strand=$strand, loc=$loc, end_pos=$end_pos;",
			"existing info: len=", $seqInfo{$seq}->[0], " strand=", $seqInfo{$seq}->[1], " loc=", $seqInfo{$seq}->[2], " end_pos=", $seqInfo{$seq}->[3];
		}
	}
	
	if(!defined $seqCnts{$seq}->{$soap}){
		$seqCnts{$seq}->{$soap} = 1;
	}else{
		$seqCnts{$seq}->{$soap} += 1;
	}
	if(!defined $seqCnts{$seq}->{"total"}){
		$seqCnts{$seq}->{"total"} = 1;
	}else{
		$seqCnts{$seq}->{"total"} += 1;
	}
	
	if($strand eq '+'){
		$num_pos++;
	}elsif($strand eq '-'){
	    $num_neg++;
	}else{
		die "strand $strand does not match pattern";
	}
}

print "In File ", $soap, "\n";
print "Number of reads mapped to luc: ", ($num_pos+$num_neg), 
      " on positive strand: $num_pos, on negative strand: $num_neg\n";
print "length distribution: \n";
foreach my $i(sort {$a<=>$b} keys %reads){
	print $i, "\t", $reads{$i}, "\n";
}
}
print "Mapped reads in each library:\n";
print join("\t", ("Seq", "Length", "Strand", "Start", "End", @ARGV)), "\n";


foreach my $seq(sort {$seqCnts{$b}->{"total"}<=>$seqCnts{$a}->{"total"}} keys %seqCnts){
	print join("\t", ($seq, @{$seqInfo{$seq}}));
	foreach my $f(@ARGV){
		if(defined $seqCnts{$seq}->{$f}){
			print "\t", $seqCnts{$seq}->{$f};
		}else{
			print "\t", 0;
		}
	}
	print "\n";
}


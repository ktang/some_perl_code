#!/usr/bin/perl -w
# from soap output, count number of small RNAs that map to luc region
use strict;

my $usage = "$0 <soap_out> [soap_out 2 ....]";
die $usage unless(@ARGV >= 1);
my ($start, $end) = (745, 2397);

foreach my $soap(@ARGV){

open(IN, $soap) or die "Can't open $soap: $!";
my %reads;
my %seqInfo;
my ($num_pos, $num_neg) = (0,0);
while(<IN>){
	my ($id, $seq, $qual, $hits, $ab, $len, $strand, $chr, $loc, $type) = split /\t/;
    my $end_pos = $loc+$len-1;
	if($loc > $end || $end_pos < $start){
		next;
	}

	$reads{$len}++;
	if(!defined $seqInfo{$seq}){
		$seqInfo{$seq} = [1, $len, $strand, $loc, $end_pos];
	}else{
		$seqInfo{$seq}->[0] += 1;
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
print "Information on top 5 most abundant reads:\n";
print join("\t", ("Seq", "Frequency", "Length", "Strand", "Start", "End")), "\n";
my $num_reads_to_show = 5;
my $n = 0;
foreach my $seq(sort {$seqInfo{$b}->[0]<=>$seqInfo{$a}->[0]} keys %seqInfo){
	$n++;
	print join("\t", ($seq, @{$seqInfo{$seq}})), "\n";
	last if($n > $num_reads_to_show);
}
}


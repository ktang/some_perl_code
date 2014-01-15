#!/usr/bin/perl -w
# check consistency between position files and brat output file
# can check position, strand only
use strict;

my $pos_file = shift @ARGV;
open(PF, $pos_file) or die "Can't open $pos_file: $!";
my %poss;
while(<PF>){
	chomp;
	next if(/pos_id/);
	my ($pos_id, $genome, $chr, $pos, $strand, $type) = split /\t/;
	$poss{$chr}->{$pos} = $strand;
}

foreach my $brat(@ARGV){
open(BF, $brat) or die "Can't open $brat: $!";
while(<BF>){
	chomp;
	my ($chr, $pos1, $pos2, $cover, $num_T, $strand) = split /\t/;
	$pos1 += 1; # brat coordinates start with 0, not 1
	if(defined $poss{$chr}->{$pos1}){
		if($poss{$chr}->{$pos1} eq $strand){
			$poss{$chr}->{$pos1} = '0';
		}else{
			print "chr $chr pos $pos1 has a different strand\n";
		}
	}else{
		print "chr $chr pos $pos1 not defined in position file\n";
	}
}
}

foreach my $chr(keys %poss){
	my %p = %{$poss{$chr}};
	foreach my $pos(keys %p){
		if($poss{$chr}->{$pos} ne '0'){
			print "chr $chr pos $pos in position file not found\n";
		}
	}
}

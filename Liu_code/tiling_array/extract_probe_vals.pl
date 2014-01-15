#!/usr/bin/perl -w
# calculate probe intensity difference
use strict;
#use Statistics::Descriptive;

my $debug = 0;

my $usage = "$0 <tpmap file> <cel files>";
die $usage unless(@ARGV >= 2);
my $tpmap = shift @ARGV;
open(TP, $tpmap) or die "Cannot open $tpmap:$!";
my %chrs;
my %probePos;
my %orderedPos;
while(<TP>){
	next if(/\#/);
	chomp;
	my @temp = split; # seq, strand, chr, start, pmx, pmy, mmx, mmy
	if(($temp[2] ne "chrM") && ($temp[2] ne "chrC")){ # not consider chloroplast or mitochondria
		 $probePos{$temp[4]}->{$temp[5]} = $temp[3];
		 push @{$orderedPos{$temp[2]}}, $temp[3];
		 $chrs{$temp[4]}->{$temp[5]} = $temp[2];
		 #if($temp[2] =~ /chr(\d)/){
		#	$chrs{$temp[4]}->{$temp[5]} = "chr" . $1;
		 #}else{
		#	die "chr ", $temp[2], " does not match pattern";
		#}
    }
}

my %values;#file->chr->start=value
my $probeLen = 25;
foreach my $cel(@ARGV){
	open(CF, $cel) or die "Can't open $cel: $!";
	if($debug){
		print STDERR "Doing cel file: $cel\n";
	}
	while(<CF>){
		chomp;
		my ($x, $y, $v) = split /\t/;
		if(defined $probePos{$x}->{$y}){
			my $chr = $chrs{$x}->{$y};
			my $start = $probePos{$x}->{$y};
			$values{$cel}->{$chr}->{$start} = log2($v);
		}else{
		#	die "probe position not defined for x=$x, y=$y";
		}
	}
}

my @files = @ARGV;
print join("\t", ("Chr", "Start", @files)), "\n";

foreach my $chr(sort keys %orderedPos){
	foreach my $p(@{$orderedPos{$chr}}){
		print join("\t", ($chr, ($p+1)));
		foreach  my $file(@files){
			print "\t", $values{$file}->{$chr}->{$p};
		}
		print "\n";
	}
}

sub log2 { 
	my $n = shift; 
	return (log($n)/log(2)); 
}

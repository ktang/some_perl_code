#!/usr/bin/perl -w
# create a flat file that contain all unique probes. Each probe has its own
# unit. The purpose is to make sure that each probe is normalized.

use strict;

my $debug = 0;
my $prog_report = 0;
my $usage = "$0 <tpmap file for unique probes>";
die $usage unless(@ARGV >= 1);
my $tpmapFile = $ARGV[0];

# print header
print join("\t", ("Probe_ID","X","Y","Probe_Sequence","Group_ID", "Unit_ID")), "\n";
my $probe_id = 1;

open(TP, $tpmapFile) or die "Cannot open $tpmapFile:$!";
my %ids;
while(<TP>){
	next if(/\#/);
	chomp;
	my @temp = split; # seq, strand, chr, start, pmx, pmy, mmx, mmy
	if(($temp[2] ne "chrM") && ($temp[2] ne "chrC")){ # not consider chloroplast or mitochondria
		my $id = $temp[2] . "_" . $temp[3];
		if(!defined $ids{$id}){
	    print join("\t", ($probe_id, $temp[4], $temp[5], $temp[0], $id, $id)), "\n";
		$ids{$id} = 1;
		$probe_id++;
		}else{
			die "id $id is seen more than once";
		}
	}
}



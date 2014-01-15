#!/usr/bin/perl -w

use strict;

my $usage = "$0 <dir> <file>";
die $usage unless(@ARGV == 2 );
my $indir = $ARGV[0];
my $infile = $ARGV[1];

#opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my $pre = "NONE";
if ($infile =~ /(\S+).soapout$/){
	$pre = $1;
}

my $output = $pre."_hit_number_profile.txt";

my $soapout = $pre."_low_hit.soap";

if (-e "$indir/$output" or $pre eq "NONE" or -e "$indir/$soapout") {
	die "$output exists or wrong name!!!\n";
}

open (IN,"$indir/$infile") or die "cannot open $infile";
open (OUT,">$indir/$output") or die "cannot open $output";
open (LOW, ">$indir/$soapout") or die "cannot open $soapout";
#six counter.
my ($one,$two, $three, $more, $one_no_mismatch, $one_mismatch)=
	(0,0,0,0,0,0);

my $last = "NONE";

while(<IN>){
	chomp;
	my @a = split "\t";
	if ($a[0] ne $last){
		$last = $a[0];
		if($a[3] == 1){
			$one++;
			print LOW $_, "\n";
			if($a[9] == 0){
				$one_no_mismatch++;
			}else{
				$one_mismatch++;
			}
		}elsif($a[3] == 2){
			print LOW $_, "\n";
			$two++;
		}elsif($a[3] == 3){
			print LOW $_, "\n";
			$three++;
		}elsif($a[3] > 3){
			$more++;
		}else{
			print STDERR $_,"\n";
			die "wrong number";
		}
	}else{
		if ($a[3] <= 3){
			print LOW $_, "\n";
		}
	}
}

print OUT join("\t",("library",">3_location","3_loacation","2_loc","1_loc" ,"perfect_match", "1_hit_with_mismatches")), "\n";
print OUT join("\t", ($pre,$more, $three, $two, $one, $one_no_mismatch, $one_mismatch)), "\n";

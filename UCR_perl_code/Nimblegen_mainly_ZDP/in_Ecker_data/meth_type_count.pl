#!/usr/bin/perl -w

use strict;

my $usage = "<$0> <input>";

die $usage unless(@ARGV == 1);


open (IN, "<$ARGV[0]")
	or die "cannot open input file:$!";
	
my $line;
my @pts;
my ($cg,$chh,$chg) =(0,0,0);
my $total = 0;
while($line = <IN>)
{
	$total++;
	chomp $line;
	@pts = split "\t", $line;
	if ($pts[3] eq "CG"){
		$cg++;
	}
	elsif ($pts[3] eq "CHG"){
		$chg++;
	}
	elsif ($pts[3] eq "CHH"){
		$chh++;
	}
	else{
		print "$line\n";
		die"unknown";
	}
}

close(IN);

print "total = $total\nCG = $cg\t",$cg/$total,"\nCHH = $chh\t",$chh/$total."\nCHG = $chg\t",$chg/$total,"\n",$cg+$chh+$chg,"\n";

print STDERR "\a";
exit;

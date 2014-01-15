#!/usr/bin/perl -w

# for the db_file 
# in (/Volumes/Macintosh_HD_2/idm1_new_met_data/UCR_server_data/Oct_1/P001)

# count the number of CG, CHG and CHH cross hyper and hypo

use strict;

my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
}

my $usage = "$0 <in_file> STDOUT";

die $usage unless (@ARGV == 1);

my $infile = shift;

die "infile" unless (-e $infile);

#my ($cg_hyper, $chg_hyper, $chh_hyper, $cg_hypo, $chg_hypo, $chh_hypo)
#	 = (0, 0, 0, 0, 0, 0);

my (%hypers, %hypos);

($hypers{CG}, $hypers{CHG}, $hypers{CHH}, $hypos{CG}, $hypos{CHG}, $hypos{CHH})
= (0, 0, 0, 0, 0, 0);
open (IN, $infile) or die "cannot open $infile: $!";
#	0		1	  2		  3		  4			5			  6		  7		  8		  9		  10	  11	  12			13
# chr1    4472    +       CG      53      0.301887        1       11      0       0       16      0       0       0.00075797

my $total = 0;


while(<IN>){
	$total++;
	chomp;
	my @a = split "\t";
	my $type = $a[3];
	my $sum = $a[7] + $a[10];
	if($sum > 0){
		if (  (   ($a[8] * $a[7] + $a[11] * $a[10]) / $sum ) > $a[5] ){
			#${$ref_hyper}{$a[0]}->{$a[1]} = 1;
			$hypers{$type}++;
		}
		else{
			$hypos{$type}++;
			#${$ref_hypo}{$a[0]}->{$a[1]} = 1;
		}
	}

	else{
		$hypos{$type} ++;
		#${$ref_hypo}{$a[0]}->{$a[1]} = 1;
	}
}

close(IN);

print STDERR $infile, ":\n";

print STDERR join("\t", ("total", "CG_hyper", "CHG_hyper", "CHH_hyper", "CG_hypo", "CHG_hypo", "CHH_hypo")), "\n";
print STDERR join("\t", ($total, $hypers{CG}, $hypers{CHG}, $hypers{CHH}, $hypos{CG}, $hypos{CHG}, $hypos{CHH} )), "\n";
print STDERR join("\t", cal_per($total, $hypers{CG}, $hypers{CHG}, $hypers{CHH}, $hypos{CG}, $hypos{CHG}, $hypos{CHH})), "\n\n";

my $sum = $hypers{CG} + $hypers{CHG} + $hypers{CHH} + $hypos{CG} + $hypos{CHG} + $hypos{CHH} ;
if($sum != $total) {
	print STDERR "\n\n\n\nwrong total\n\n\n\n";
}

exit;

sub cal_per{
	my ($t) = shift @_;
	my @a   = @_;
	my @b   = ();
	for my $i (0..5){
		$b[$i] = sprintf ("%.1f", 100 * $a[$i] / $t) . "%";
	}
	return ("100%", @b);
}



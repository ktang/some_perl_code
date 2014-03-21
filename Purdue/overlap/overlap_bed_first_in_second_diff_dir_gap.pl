#!/usr/bin/perl -w
# two bed files as input
# foreach region in first file, check whether there is a region in overlap with it
# output
use strict;

my $debug = 0;
if ($debug){
	print STDERR "debug: $debug\n";
}

print STDERR "input must have a head with Start in it.\n\n";

my $usage = "$0 <bed_file1> <bed_file2> <gap>";

die $usage unless (@ARGV == 3);

#my ($indir, $input1, $input2, $allowed_gap) = @ARGV[0..3];

my $in1  = shift or die;
my $in2  = shift or die;
#my $allowed_gap = 0;
my $allowed_gap = shift;

die unless (defined $allowed_gap);

die "wrong directory\n\n" unless (-e $in1 and -e $in2);

open (IN1, $in1) or die "cannot open $in1: $!";
open (IN2, $in2) or die "cannot open $in2: $!";
#open (OUT, ">$output") or die "cannot open $output: $!";

my %beds;
my $num_f2 = 0;
my $h = <IN2>;
while (<IN2>){
	next if (/^#/);
	$num_f2++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $new_start = $start - $allowed_gap;
	if($new_start <=0 ){$new_start = 1 }
	
	for my $i ( $new_start..($end + $allowed_gap)){
		$beds{$chr}->[$i] = 1;
	} 
}

my $num_f1 = 0;
my $num_overlap = 0;
$h = <IN1>;
while (<IN1>){
	next if (/^#/);
	$num_f1 ++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $flag = 0;
	for my $i ($start..$end){
		if(defined $beds{$chr}->[$i]){
			$flag = 1;
			last;
		}
	}
	
	if ($flag == 1){$num_overlap ++}
	
}

my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);

print "$in2: $num_f2\n";
print "$in1: $num_f1\n";
print "overlap (gap:$allowed_gap): ",  "$num_overlap / $num_f1 = $per%\n\n";
exit;
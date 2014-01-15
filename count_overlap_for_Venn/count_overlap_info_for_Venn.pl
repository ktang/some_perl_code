#!/usr/bin/perl -w
use strict;
use File::Spec;
my $debug = 0;
my $usage = "$0 \n <input>  <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
#my $sample_index = shift or die; 
die unless (-e $input);

my $output = shift or die;
die if (-e $output);

open(IN, $input) or die;
my $h = <IN>;
chomp $h;
my @h_a = split "\t", $h;

my @labels;
my @indexs;
foreach my $i(0..$#h_a){
	if($h_a[$i] =~ /^overlap_(\S+)/){
		push @labels, $1;
		push @indexs, $i;
	}
}

if ($debug) {
	foreach my $i(0..$#indexs){
		print STDERR join("\t", ($indexs[$i], $labels[$i])),"\n";
	}
	exit;
}

open(OUT, ">>$output") or die ;


my $total = 0;
#my $power = 2**$sample_index;
#my @numbers = (0) x $power ;
my %nums;
while (<IN>){
	$total++;
	chomp;
	my @a = split "\t";
	my @flags = ();
	for my $i (0..$#indexs){
		#my $base = 2 ** ($i - 1);
		if($a[$indexs[$i]] =~/chr\d_\d/){
			push @flags, "Overlap";
		}else{
			die $_ unless ( $a[$indexs[$i]]  eq "NOT" );
			push @flags, "Not";
		}
	}
	my $flag = join("\t", @flags);
	$nums{$flag}++;
}

close IN;

my $sum = 0 ;
print OUT join("\t", ("num", @labels)), "\n";
foreach my $k (sort keys %nums){
	$sum += $nums{$k};
	print OUT join("\t", ($nums{$k}, $k)), "\n";
}

die "sum != total, sum= $sum, total = $total" if ($sum != $total);

exit;
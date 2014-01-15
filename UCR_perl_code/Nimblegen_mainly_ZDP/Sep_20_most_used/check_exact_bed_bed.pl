#!/usr/bin/perl -w

use strict;

#short for bed, long for gff
my $usage = "$0 <bed_file1> <bed_file2> ";
die $usage unless(@ARGV == 2);

my $fname1=$ARGV[0]; #input filenam
my $fname2=$ARGV[1]; #second filename

my @sps;
my @lps;
my $line;
my ($send,$sstart,$lend,$lstart);
my ($long_num,$short_num,$overlap) =(0,0,0);

open (IN1,"<$fname1")
	or die "cannot open file $fname1:$!";
open (IN2,"<$fname2")
	or die "cannot opne file $fname2:$!";

my @beds1 = <IN1>;
close(IN1);

my @beds2 = <IN2>;
close (IN2);

my $num_same = 0;

LOOP: foreach my $peak1 (@beds1)
{
	my $findflag = 0;
	chomp $peak1;
	my @pts1 = split "\t", $peak1;
	foreach my $peak2 (@beds2)
	{
		chomp $peak2;
		my @pts2 = split "\t", $peak2;
		if ( ($pts1[0] eq $pts2[0]) and  ($pts1[1] eq $pts2[1]) and ($pts1[2] eq $pts2[2]))
		{
			$findflag = 1;
			$num_same++;
			next LOOP;
		} 
	}	
	if ($findflag == 0)
	{
		print STDERR "$peak1\n";
	}
}

print STDERR "$fname1 contains $#beds1+1 line\n$fname2 contains $#beds2+1 line\n";
print STDERR "the same is $num_same\n";

print STDERR "\a";
exit;
#!/usr/bin/perl -w

# calculate the mean of the two sorted sgr files

use strict;

my $usage = "$0 <input1> <input2>";

die $usage unless(@ARGV == 2);

open (IN1 , "<$ARGV[0]")
  or die "Cannot open input_file1: $!";
  
open (IN2,  "<$ARGV[1]")
  or die "Cannot open input_File2: $!";
  
  
  
  my $l1;
  my $l2;
  my @p1;
  my @p2;
  my $mean;
  
#  print OUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup\n";
 # $l1 = <IN1> ,$l2 = <IN2>;
  
  while($l1 = <IN1> ,$l2 = <IN2>)
  {
	chomp $l1;
	chomp $l2;
	@p1 = split /\t/, $l1;
	@p2 = split /\t/, $l2;
	if (($p1[0] eq $p2[0]) and ($p1[1] eq $p2[1]))
	{
		$mean = 0.5*($p1[2] + $p2[2]);
	}
	else
	{
		print STDERR "$l1\t$l2\n";
		die "wrong probe";	
	}
	$p1[2] =$mean;
	my $out_line = join "\t",@p1;
	print "$out_line\n";
	
  }

  close (IN1);
  close (IN2);
  
  exit;

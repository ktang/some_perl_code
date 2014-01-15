#!/usr/bin/perl -w

# calculate the mean of the two sorted gff files

use strict;

my $usage = "$0 <input1> <input2> <output>";

die $usage unless(@ARGV == 3);

open (IN1 , "<$ARGV[0]")
  or die "Cannot open input_file1: $!";
  
open (IN2,  "<$ARGV[1]")
  or die "Cannot open input_File2: $!";
  
open (OUT,">$ARGV[2]")
  or die "Cannot open output files: $!";
  
  
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
	$mean = 0.5*($p1[5] + $p2[5]);
	print OUT "$p1[0]\t$p1[1]\t$p1[2]\t$p1[3]\t$p1[4]\t$mean\t$p1[6]\t$p1[7]\t$p1[8]\n";
	
  }

  close (IN1);
  close (IN2);
  close (OUT);
    
  print "\a";
  exit;
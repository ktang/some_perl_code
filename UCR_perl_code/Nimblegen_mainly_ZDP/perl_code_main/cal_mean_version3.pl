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
	if ($p1[8] eq $p2[8])
	{
		$mean = 0.5*($p1[5] + $p2[5]);
	}
	else
	{
		die "wrong probe";	
	}
	$p1[5] =$mean;
	my $out_line = join "\t",@p1;
	print OUT "$out_line\n";
	
  }

  close (IN1);
  close (IN2);
  close (OUT);
    
  print "\a";
  exit;
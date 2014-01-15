#!/usr/bin/perl -w

# calculate the difference between mutants and WT
# note only Col0 #2 is used

use strict;

my $usage = "$0 <C2> <input> <output>";

die $usage unless(@ARGV == 3);

open (WT , "<$ARGV[0]")
  or die "Cannot open wt: $!";
  
open (IN,  "<$ARGV[1]")
  or die "Cannot open input_File: $!";
  
open (OUT,">$ARGV[2]")
  or die "Cannot open output files: $!";
  
  my $l1;
  my $l2;
  my @p1;
  my @p2;
  my $diff;

  while($l1 = <IN> ,$l2 = <WT>)
  {
	chomp $l1;
	chomp $l2;
	@p1 = split /\t/, $l1;
	@p2 = split /\t/, $l2;
	if ($p1[8] eq $p2[8])
	{
		$diff = $p1[5]-$p2[5];
	}
	else 
	{
		die "not the same";
	}
	$p1[5] = $diff;
	my $out_line = join "\t",@p1;
	print OUT "$out_line\n";
  }

  close (IN);
  close (WT);
  close (OUT);
    
  print "\a";
  exit;
#!/usr/bin/perl -w

# calculate the difference between mutants and WT
# note only Col0 #2 is used

use strict;

my $usage = "$0 <col0> <mutant> <output>";

die $usage unless(@ARGV == 3);

open (WT , "<$ARGV[0]")
  or die "Cannot open wt: $!";
  
open (MU,  "<$ARGV[1]")
  or die "Cannot open input_File: $!";
  
open (OUT,">$ARGV[2]")
  or die "Cannot open output files: $!";
  
  my $l1;
  my $l2;
  my @p1;
  my @p2;
  my $diff;

  while($l1 = <MU> ,$l2 = <WT>)
  {
	chomp $l1;
	chomp $l2;
	@p1 = split /\t/, $l1;
	@p2 = split /\t/, $l2;
	if (($p1[0] eq $p2[0]) and ($p1[3] eq $p2[3]))
	{
		$diff = $p1[7]-$p2[7];
	}
	else 
	{
		die "not the same";
	}
	$p1[5] = $diff;
	$p1[6] = ".";
	$p1[7] = ".";
	my $out_line = join "\t",@p1;
	print OUT "$out_line\n";
  }

  close (MU);
  close (WT);
  close (OUT);
    
  print STDERR "\a";
  exit;
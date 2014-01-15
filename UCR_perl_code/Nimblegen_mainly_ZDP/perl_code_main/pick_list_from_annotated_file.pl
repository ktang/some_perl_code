#!/usr/bin/perl -w

=head2
not finished
this script is used as pick trimmeaned part of
list from the original file
=cut

use strict;

my $usage = "$0 <input1> <input2> <output>";

die $usage unless(@ARGV == 3);

open (BED , "<$ARGV[0]")
  or die "Cannot open input_file1: $!";
  
open (ANNO,  "<$ARGV[1]")
  or die "Cannot open input_File2: $!";
  
open (OUT,">$ARGV[2]")
  or die "Cannot open output files: $!";
  
  
  my $l1;
  my $l2;
  my @p1;
  my @p2;
  my $mean;
  my %regions;

  
  while($l1 = <BED>)
  {
	chomp $l1;
	@p1 = split /\t/, $l1;
	my $chr = lc($p1[0]);
	$regions{$p1[0]}->{$p1[1]} = 1;
	
  }
  
  while

  close (BED);
  close (ANNO);
  close (OUT);
    
  exit;
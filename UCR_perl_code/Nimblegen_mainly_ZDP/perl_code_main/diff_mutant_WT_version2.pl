#!/usr/bin/perl -w

# calculate the difference between mutants and WT

use strict;

my $usage = "$0 <input>";

die $usage unless(@ARGV == 1);

open (IN , "<$ARGV[0]")
  or die "Cannot open input_file1: $!";
  
open (WT,  "<Col0_mean.gff")
  or die "Cannot open input_File2: $!";
  
   
  $ARGV[0] =~ /(\S+)_mean\.gff/;
  my $pre = $1;
  my $output = $1."_diff_form_WT.gff";

  
open (OUT,">$output")
  or die "Cannot open output files: $!";
  
  
  my $l1;
  my $l2;
  my @p1;
  my @p2;
  my $diff;
  
#  print OUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup\n";
 # $l1 = <IN> ,$l2 = <WT>;
  
  while($l1 = <IN> ,$l2 = <WT>)
  {
	chomp $l1;
	chomp $l2;
	@p1 = split /\t/, $l1;
	@p2 = split /\t/, $l2;
	$diff = $p1[5]-$p2[5];
	print OUT "$p1[0]\t$p1[1]\t$p1[2]\t$p1[3]\t$p1[4]\t$diff\t$p1[6]\t$p1[7]\t$p1[8]\n";
	
  }

  close (IN);
  close (WT);
  close (OUT);
    
  print "\a";
  exit;
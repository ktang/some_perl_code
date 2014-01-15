#!/usr/bin/perl -w

# use sorted_exact_probe_position file to sort
# the gff file

use strict;

my $usage = "$0 <input_file>";

die $usage unless(@ARGV == 1);

open (SORTED,  "</Users/kaitang/Desktop/Nimblegen/real_data/probe_ID_position/sorted_exact_probe_position")
  or die "Cannot open sorted_File: $!";
  
open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";
  
  
  $ARGV[0] =~ /(\S+)\.gff/;
  my $pre = $1;
  my $output = $1."_sorted.gff";
open (OUT,  ">$output")
  or die "Cannot open output_File: $!";

 
  my $line;
  my @parts;
  my %value;
  
  while($line = <IN>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	
	$value{$parts[2]} = $parts[5];
  }

  close (IN);
  
  
  while($line =<SORTED>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	print OUT $line,"\t$value{$parts[2]}\t.\t.\t.\n";
	
  }
  
  close(SORTED);
  close(OUT);
  
  print "\a";
  exit;
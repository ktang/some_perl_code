#!/usr/bin/perl -w

# use sorted_probe_position file, filter the oringinal 
#file

use strict;

my $usage = "$0 <input_file>";

die $usage unless(@ARGV == 1);

open (SORTED,  "/Users/kaitang/Desktop/Nimblegen/probe_files/sorted_position_ID_v2")
  or die "Cannot open sorted_File: $!";
  
open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";
  
  
  $ARGV[0] =~ /^(.+)\.txt/;
  my $pre = $1;
  my $output = $1."_filter.txt";
  
  
open (OUT,  ">$output")
  or die "Cannot open output_File1: $!";

 
  my $line;
  my @parts;
  my %v1;
  my %v2;
  
  my $header = <IN>;
  
  while($line = <IN>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	
	$v1{$parts[0]} = $parts[1];
    $v2{$parts[0]} = $parts[2];
	
	
  }

  close (IN);
  
  print OUT $header;
  
  while($line =<SORTED>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	if (exists $v1{$parts[3]})
	{
	print OUT "$parts[3]\t$v1{$parts[3]}\t$v2{$parts[3]}\n";
		}
  }
  
  close(SORTED);
  close(OUT);
  
  print "\a";
  exit;
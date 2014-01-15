#!/usr/bin/perl -w

# use sorted_exact_probe_position file, convert the normalized value
# writed out from the R to the file used for MATALg

use strict;

my $usage = "$0 <input_file>";

die $usage unless(@ARGV == 1);

open (SORTED,  "./sorted_exact_probe_position")
  or die "Cannot open sorted_File: $!";
  
open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";
  
  
  $ARGV[0] =~ /^(...)/;
  my $pre = $1;
  my $output1 = $1."_Gquantile_sample1_sorted.gff";
  my $output2 = $1."_Gquantile_sample2_sorted.gff";
  
open (OUT1,  ">$output1")
  or die "Cannot open output_File1: $!";

open (OUT2,  ">$output2")
  or die "Cannot open output_File2: $!"; 
 
  my $line;
  my @parts;
  my %v1;
  my %v2;
  
  <IN>;
  
  while($line = <IN>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	
	$v1{$parts[0]} = $parts[1];
    $v2{$parts[0]} = $parts[2];
	
	
  }

  close (IN);
  
  
  while($line =<SORTED>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	if (exists $v1{$parts[2]})
	{
	print OUT1 $line,"\t$v1{$parts[2]}\t.\t.\t.\n";
	print OUT2 $line,"\t$v2{$parts[2]}\t.\t.\t.\n";
	}
  }
  
  close(SORTED);
  close(OUT1);
  close(OUT2);
  print "\a";
  exit;
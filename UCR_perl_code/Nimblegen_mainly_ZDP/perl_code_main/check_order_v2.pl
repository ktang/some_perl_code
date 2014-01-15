#!/usr/bin/perl -w

# check whether probe is ordered
# by the position in chr

use strict;

my $usage = "$0 <input_gff>";

die $usage unless(@ARGV == 1);

open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";
  
  
  my @lines;
  
  @lines = <IN>;
  close(IN);

 
  my @ps1;
  my @ps2;
  my $count = 0;
  my $i = 0;
  
  for($i = 0; $i < $#lines; $i++)
  {
  		@ps1 = split "\t", $lines[$i];
  		@ps2 = split "\t", $lines[$i+1];
  		if( ($ps1[0] eq $ps2[0] )  and ( $ps1[1] > $ps2[1] ))
  		{
  			#print "$lines[$i]\n$lines[$i+1]";
  			$count++;
  			}
  }
    
    
  print "$count\n";
  print "\a";
  exit;
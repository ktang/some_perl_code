#!/usr/bin/perl -w

# check whether probe is ordered
# by the position in chr

use strict;

my $usage = "$0 <input_gff>";

die $usage unless(@ARGV == 1);

open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";

  my %position;
  my $line;
  my @parts;
  my $former = 0;
  my $former_probe;
    
  $line = <IN>;
  while($line = <IN>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	
	if($parts[3] <= $former)
	{
		print $line,"\t$former\t$former_probe\n";
		$former = $parts[3];	
	}	
	else
	{
		$former = $parts[3];
	}
	$former_probe = $parts[2];
  }

  close (IN);
  
  
  print "\a";
  exit;
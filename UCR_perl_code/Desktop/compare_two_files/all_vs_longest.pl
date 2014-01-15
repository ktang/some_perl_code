#!/usr/bin/perl -w

# get the reads which mapped to 
#the all_transcript but not the longest


use strict;

my $usage = "$0 <longest> <all> <output>";

die $usage unless(@ARGV == 3);

open (LONG,  "<$ARGV[0]")
  or die "Cannot open $ARGV[0]: $!";

open (ALL, "<$ARGV[1]")
  or die "Cannot open $ARGV[1]: $!";
  
  open (OUT, ">$ARGV[2]")
  	 or die "Cannot open $ARGV[2]: $!";

  
  my $line;
  my @parts;
  my %reads;
  
  while($line = <LONG>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	$reads{$parts[0]} = -1;
  }

  close (LONG);
  
  while ($line = <ALL>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	if( exists $reads{$parts[0]})
	{}
	else 
	{
		print OUT "$line\n";	
	}
  }
  
 
  close (ALL);
    close (OUT);
  

  print "\a";
  exit;
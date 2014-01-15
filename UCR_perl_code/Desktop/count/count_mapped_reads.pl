#!/usr/bin/perl -w

# input the soap output (which the seq is not uniq with RDN number)
# and count how many reads has been mapped to the ref seq

my $usage = "$0 <soapout>";
die $usage unless(@ARGV == 1);

open (IN, "<$ARGV[0]")
  or die "Cannot open inFile: $!";


  my @parts;
  my %reads;
  my $count = 0;
  my $line;
  
  while($line = <IN>)
  {
     
    @parts = split /\t/,$line;
    if (exists $reads{$parts[0]})
    {
      ;
    }
    else
    {
    	$reads{$parts[0]} = -1;
    	$count++;
    }
    
  }
    close (IN);
    print "\ntotal mapped seq is:$count\n";
    print "\a";
  exit;
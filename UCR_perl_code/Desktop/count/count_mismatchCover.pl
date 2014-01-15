#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2>";
die $usage unless(@ARGV >= 2);
my ($inFile1, $inFile2, $outFil) = @ARGV[0..1];

my $IN1;
my $IN2;

open ($IN1, '<', $inFile1)
  or die "Cannot open $inFile1: $!";
   
open ($IN2, '<', $inFile2)
  or die "Cannot open $inFile2: $!";


  

  my %uniq;
  my $line;
  my @parts;
  my $count = 0;
  
  while($line = <$IN1>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$uniq{$parts[0]} = -1;
  }
  
  while($line = <$IN2>)
  {
     chomp($line);
     @parts = split /\t/,$line;
   if (exists $uniq{$parts[0]})
    {
      next;
    }
    else
    {
      
      $count += $parts[1];
    } 
  }



  close ($IN1);
  close ($IN2);
  print "total number of reads is:\t",$count,"\n";
  print "\a";
  exit;
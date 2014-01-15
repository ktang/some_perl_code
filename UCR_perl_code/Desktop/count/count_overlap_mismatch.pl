#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2> <in3> <in4>";
die $usage unless(@ARGV >= 4);
my ($inFile1, $inFile2, $inFile3, $inFile4) = @ARGV[0..3];

my $IN1;
my $IN2;
my $IN3;
my $IN4;

open ($IN1, '<', $inFile1)
  or die "Cannot open $inFile1: $!";
   
open ($IN2, '<', $inFile2)
  or die "Cannot open $inFile2: $!";

open ($IN3, '<', $inFile3)
  or die "Cannot open $inFile3: $!";

open ($IN4, '<', $inFile4)
  or die "Cannot open $inFile4: $!";

  

  my %uniq1;
  my %uniq2;
  my %uniq3;
  my $line;
  my @parts;
  my $count = 0;
  
  while($line = <$IN1>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$uniq1{$parts[0]} = -1;
  }
  
   while($line = <$IN2>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$uniq2{$parts[0]} = -1;
  }
  
  
   while($line = <$IN3>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$uniq3{$parts[0]} = -1;
  }
  
  
  
  while($line = <$IN4>)
  {
     chomp($line);
     @parts = split /\t/,$line;
   if ( (exists $uniq3{$parts[0]})  && (! (exists $uniq1{$parts[0]})) && (! (exists $uniq2{$parts[0]}))  )
    {
      $count += $parts[1];
    }
  }



  close ($IN1);
  close ($IN2);
  close ($IN3);
  close ($IN4);
  print "total number of overlap is:\t",$count,"\n";
  print "\a";
  exit;
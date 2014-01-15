#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2>";
die $usage unless(@ARGV >= 2);
my ($inFile1, $inFile2, $outFile) = @ARGV[0..2];

my $IN1;
my $IN2;
my $OUT;

open ($IN1, '<', $inFile1)
  or die "Cannot open $inFile1: $!";
   
open ($IN2, '<', $inFile2)
  or die "Cannot open $inFile2: $!";
  
 open($OUT, '>', $outFile)
	or die "Cannot open $outFile: $!"; 


  

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
      
      print $OUT $line,"\n";
    } 
  }



  close ($IN1);
  close ($IN2);
  close ($OUT);
  
  print "\a";
  exit;
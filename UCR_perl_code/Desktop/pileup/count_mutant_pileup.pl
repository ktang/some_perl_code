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
 
  
  while($line = <$IN1>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$uniq{$parts[0].$parts[1]} = $parts[5];
  }
  
  while($line = <$IN2>)
  {
     chomp($line);
     @parts = split /\t/,$line;
   if (exists $uniq{$parts[0].$parts[1]})
    {
       if( ($parts[5] > 0.5&& $uniq{$parts[0].$parts[1]} < 0.5)
       || ($parts[5] < 0.5&& $uniq{$parts[0].$parts[1]} > 0.5))
     {
       print  $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t",$parts[5],"\t",$uniq{$parts[0].$parts[1]},
		"\n";
      
      }
    }
  
    
  }



  close ($IN1);
  close ($IN2);
  
  print "\a";
  exit;
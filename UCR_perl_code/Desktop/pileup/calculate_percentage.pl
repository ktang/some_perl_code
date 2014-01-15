#!/usr/bin/perl -w

use strict;

my $usage = "$0 <in> <out>";
die $usage unless(@ARGV >= 2);
my ($inFile, $outFile) = @ARGV[0..1];

my $IN;
my $OUT;

open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";

open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  
  my $i;
  my $count;
  my $line;
  my @parts;
  my $char;
  my $len;
  my ($num_A, $num_T, $num_C, $num_G);
  while($line = <$IN>)
  {
	chomp($line);
	$count = 0;
	($num_A, $num_T, $num_C, $num_G) =(0,0,0,0);
	@parts = split /\t/,$line;

=POD
      if ($parts[4] =~ /\+/ || $parts[4] =~ /-/)
      {
      print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t","-1","\n";
      	next;
      }
=cut

      $len = length($parts[4]);
      
    for ($i = 0 ; $i < $len; $i++)
      {
      	$char = substr($parts[4], $i,1 );
   		if ( $char eq "\." ||  $char eq ",")
      	{
      		$count++;
      	}
      	elsif($char eq "a" || $char eq "A")
      	{
      		$num_A++;
      	}
      	elsif ($char eq "t" || $char eq "T")
      	{
      		$num_T++;
      	}
      	elsif ($char eq "c" || $char eq "C")
      	{
      		$num_C++;
      	}
      	elsif ($char eq "g" || $char eq "G")
      	{
      		$num_G++;
      	}
      }
      
      my $total = $count/$parts[3]+$num_A/$parts[3]+$num_T/$parts[3]+$num_C/$parts[3]+$num_G/$parts[3];
		print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t",$count/$parts[3],
		"\t",$num_A/$parts[3],"\t",$num_T/$parts[3],"\t",
		$num_C/$parts[3],"\t",$num_G/$parts[3],"\t",$total,"\n";
      
    }
    
   

  close ($IN);
  close ($OUT);
  print "\a";
  exit;
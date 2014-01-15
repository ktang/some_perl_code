#!/usr/bin/perl -w
# after extract_onestep.pl
#use new_find_max.pl to get the 
# subset of position with the max_mismatch


use strict;

my $usage = "$0 <pileup_extracted> <out>";
die $usage unless(@ARGV == 1);


open (IN, "<$ARGV[0]")
  or die "Cannot open inFile: $!";

open (OUT, ">$ARGV[0].max_mismatch")
  or die "Cannot open outFile: $!";
  
  my $i;
  my $total;
  my $line;
  my @parts;
  my $char;
  my $len;
  my @bp =("T","C","G");
  my $j;
  
  while($line = <IN>)
  {
	chomp($line);
	
	
	@parts = split /\t/,$line;

     
      my ($max,$max_bp) =($parts[6],"A");
      
      for ($i = 7; $i <= 9 ; $i++)
      {
      	if ($parts[$i] > $max)
      	{
      		($max,$max_bp) =($parts[$i],$bp[$i-7]);
      	}
      }
      
      $j = 0;
      for ($i = 6; $i <= 9 ; $i++)
      {
      	if ($parts[$i] == $max)
      	{
      		$j++;
      	}
      }
   	  
   	  
   	  $total = $parts[5]+ $max;
     
     if ( $j == 1 )
     {
       
      print OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
	  "\t",$parts[3],"\t",$parts[4],"\t",$parts[5],"\t",
	  $max,"\t",$max_bp,"\t",$total,"\n";
     }
     
=pod 
     elsif ($j > 1)
     {
     print parts[0],"\t",$parts[1],"\t",$parts[2],
	  "\t",$parts[3],"\t",$parts[4],"\t",$parts[5],"\t",
	  $max,"\t",$max_bp,"\t",$total,"\n";
	  }
=cut
    }
    
   

  close (IN);
  close (OUT);
  print "\a";
  exit;
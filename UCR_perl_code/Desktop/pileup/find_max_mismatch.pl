#!/usr/bin/perl -w
# batch convert fastq file to fasta file
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
  
  my $line;
  my @parts;
  my $char;
  my $len;
  my @bp =("T","C","G");
  
  while($line = <$IN>)
  {
	chomp($line);
	
	my $flag = 0;
	@parts = split /\t/,$line;

=POD
      if ($parts[4] =~ /\+/ || $parts[4] =~ /-/)
      {
      print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t","-1","\n";
      	next;
      }
       $len = length($parts[4]);
       for ($i = 0 ; $i < $len; $i++)
      {
      	$char = substr($parts[4], $i,1 );
   		if ( $char eq "\." ||  $char eq ",")
      	{
      		$count++;
      	}
      	
      }
=cut

     
      my ($max,$max_bp) =($parts[6],"A");
      for ($i = 7; $i <= 9 ; $i++)
      {
      	if ($parts[$i] > $max)
      	{
      		($max,$max_bp) =($parts[$i],$bp[$i-7]);
      	}
      	elsif ($parts[$i] == $max)
      	{
      		$flag = 1;
      	}
   	  }
     
     if ( $flag == 0 )
     {
      my $total = $parts[5]+ $max;
      print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
	  "\t",$parts[3],"\t",$parts[4],"\t",$parts[5],"\t",
	  $max,"\t",$max_bp,"\t",$total,"\n";
     }
    }
    
   

  close ($IN);
  close ($OUT);
  print "\a";
  exit;
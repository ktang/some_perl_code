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
  my $count;
  my $line;
  my @parts;
  my $char;
  my $len;
  
  while($line = <$IN>)
  {
	chomp($line);
	$count = 0;
	
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
      	
      }
      
      my $ratio = $count/$parts[3];
      
      if ( $ratio < 0.95)
      {
		print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t",$ratio,"\n";
      }
    }
    
   

  close ($IN);
  close ($OUT);
  print "\a";
  exit;
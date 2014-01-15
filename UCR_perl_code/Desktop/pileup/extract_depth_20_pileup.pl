#!/usr/bin/perl -w
# batch convert fastq file to fasta file
use strict;

my $usage = "$0 <in> <out>";
die $usage unless(@ARGV >= 2);
my ($inFile,$outFile) = @ARGV[0..1];

my $IN;
my $OUT;

open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";

open ($OUT, '>', $outFile)
	or die "Cannot open $outFile:$!";
	
  
  my $position_count = 0;
  my $line;
  my @parts;
 
  
  while($line = <$IN>)
  {
	chomp($line);
	
	@parts = split /\t/,$line;

    if ($parts[3] >= 20)
    {
     	 $position_count++;
     	 print $OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\n";
     	
     }
  }

  close ($IN);
  close ($OUT);
  print "position which depth >= 20:$position_count\n",
 
  
  print "\a";
  exit;
  
  
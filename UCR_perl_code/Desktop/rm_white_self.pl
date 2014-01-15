#!/usr/bin/perl -w
# batch convert fastq file to fasta file
use strict;

my $usage = "$0 <in>";
die $usage unless(@ARGV >= 1);
my $inFile = $ARGV[0];

my $IN;
my $OUT;

open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";

open ($OUT, '>', $inFile)
  or die "Cannot open $inFile: $!";
  my $i = 0;
  my $line;
  while($line = <$IN>)
  {


    if ($line =~ m/\sR/)
    {
      $line =~ s/(\s)R/+R/ ;
    }

    print $OUT $line;
  }

  close ($IN);
  close ($OUT);
  print "\a";
  exit;
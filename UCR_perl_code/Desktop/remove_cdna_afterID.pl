#!/usr/bin/perl -w
# batch convert fastq file to fasta file


my $usage = "$0 <in> <out>";
die $usage unless(@ARGV >= 2);
my ($inFile, $outFile) = @ARGV[0..1];

my $IN;
my $OUT;

open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";

open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  

  my $line;
  
  while($line = <$IN>)
  {


    if ($line =~ m/\.\d+(\ \|.*FORWARD)/)
    {
    
      $line =~ s/$1/why/ ;
      print $1,"\n";
    
    }

    print $OUT $line;
  }

  close ($IN);
  close ($OUT);
  print "\a";
  exit;

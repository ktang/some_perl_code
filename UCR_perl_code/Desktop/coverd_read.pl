#!/usr/bin/perl -w
# batch convert fastq file to fasta file


my $usage = "$0 <in>";
die $usage unless(@ARGV >= 1);
my $inFile = $ARGV[0];

my $IN;


open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";


  my $i;
  my $total = 0;
  my $line;
  my @parts;
  
  while($line = <$IN>)
  {

     chomp($line);
     @parts = split /\t/,$line;

     $total += $parts[1];

    
  }

  close ($IN);
  print "total number of reads is:\t",$total,"\n";
  print "\a";
  exit;

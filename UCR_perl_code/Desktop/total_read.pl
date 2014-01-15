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
  
  while($line = <$IN>)
  {


    if ($line =~ m/RDN=(\d+)/)
    {
      $i = $1 ;
      $total += $i
    }

    
  }

  close ($IN);
  print "total number of reads is:\t",$total,"\n";
  print "\a";
  exit;

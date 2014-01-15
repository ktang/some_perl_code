#!/usr/bin/perl -w
my $usage = "$0 <in> <out>";
die $usage unless(@ARGV >= 2);
my ($inFile, $outFile) = @ARGV[0..1];

my $IN;
my $OUT;

open ($IN, '<', $inFile)
  or die "Cannot open $inFile: $!";

open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  my %count;
  my %uniq;
  my $line;
  my @parts
  while($line = <$IN>)
  {
     
    @parts = split /\t/,$line;
    if (exists $uniq{$parts[8].$parts[0]})
    {
      next;
    }
    else
    {
      $uniq{$parts[8].$parts[0]} = $parts[1];
      $count{$parts[8]} += $parts[1];
    }
    
  }
  foreach my $gene(sort keys %count)
    print $OUT "$gene\t$count{$gene}\n";

  close ($IN);
  close ($OUT);
  print "\a";
  exit;
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
  my @parts;
  
  while($line = <$IN>)
  {
     
    @parts = split /\t/,$line;
    
=pod    
    if (defined $uniq{$parts[8]}->{$parts[0]})
    {print $parts[0]," map more than one position in ",$parts[8],"\n";}
    else
    {
      $uniq{$parts[8]}->{$parts[0]} = -1;
    }
=cut

    $count{$parts[8]} += $parts[1];
    
  }
  foreach my $name(sort keys %count)
  {
    print $OUT $name,"\t",$count{$name},"\n";
  }
  close ($IN);
  close ($OUT);
  print "\a";
  exit;
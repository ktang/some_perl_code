#!/usr/bin/perl -w
my $usage = "$0 <in>";
die $usage unless(@ARGV == 1);
my ($inFile) = $ARGV[0];



open (IN, '<', $inFile)
  or die "Cannot open $inFile: $!";


  my $count;
  my %uniq;
  my $line;
  my @parts;
  <IN>;
  while($line = <IN>)
  {
     
    @parts = split /\t/,$line;
    if (exists $uniq{$parts[0]})
    {
      next;
    }
    else
    {
      $uniq{$parts[0]} = -1;
      $count++;
    }
    
  }
    close (IN);
    
    print "uniq lines:","$count\n";
    print "\a";
  exit;
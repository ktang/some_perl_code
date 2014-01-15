#! /usr/bin/perl -w
#read from the TAIR cdna file, get
#the information which strand the transcript is

use strict;

my $usage = "$0 <in> <out>";
die $usage unless (@ARGV == 2);

open (IN, "<$ARGV[0]")
  or die "Cannot open $ARGV[0]: $!";

open (OUT, ">$ARGV[1]")
  or die "Cannot open $ARGV[1]: $!";

my $line;

while($line = <IN>)
{
   chomp $line;
 	if ($line =~ /^>/)
 	{
 		if( $line =~ m/^>(AT.G.......).+FORWARD$/	)
 		{
 			print OUT "$1\t+\n";
 		}
 		elsif ($line =~ m/^>(AT.G.......).+REVERSE$/)
 		{
 			print OUT "$1\t-\n";
  		}
  		
  		else
  		{
  			print $line,"\n";	
  		}
 	}
 	else
 	{;}  
}

close (IN);
close (OUT);
print "\a";
exit;
#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 17 Nov 2010

change small letter "chr" to big letter "Chr"
Nov 17,2010
=cut

use strict;
my $usage = "$0 <input_bam>";
die $usage unless(@ARGV == 1);

open (IN, "$ARGV[0]")
  or die "Cannot open input_File: $!";

my $line;

  while($line = <IN>)
  {


    if ($line =~ m/chr/)
    {
      $line =~ s/chr/Chr/g ;
    }

    print  $line;
  }

  close (IN);
 
  
  exit;

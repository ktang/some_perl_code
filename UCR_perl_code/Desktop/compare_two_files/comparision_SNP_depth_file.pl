#!/usr/bin/perl -w

# compare whether the two SNP_depth files are 
# the same. use for redo analysis

use strict;

my $usage = "$0 <former> <redo>";

die $usage unless(@ARGV == 2);

open (FORMER,  "<$ARGV[0]")
  or die "Cannot open former_File: $!";

open (RE_DO, "<$ARGV[1]")
  or die "Cannot open redo_File: $!";
  
  my %a;
  my $line;
  my @parts;
  
  while($line = <FORMER>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	$a{$parts[0].$parts[1]} = $parts[2];
  }

  close (FORMER);
  
  while ($line = <RE_DO>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	if (exists $a{$parts[0].$parts[1]} )
	{
			if($a{$parts[0].$parts[1]} == $parts[2])
			{}
			else 
			{
				print "wrong:new\t$line\told:$a{$parts[0].$parts[1]}\n";
			}
	}
	else
	{
		print "don't exists $line in old file";
	}
  }
  close (RE_DO);
  print "\a";
  exit;
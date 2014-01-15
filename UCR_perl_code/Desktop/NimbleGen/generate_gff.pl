#!/usr/bin/perl -w

# use the Arab_thal_H1_all_meth.pos file and the ten preprocessed
# files to generate the .gff files for MEDME read in.


use strict;

my $usage = "$0 <id_pos_file> <preprocessed_file> <gff_file>";

die $usage unless(@ARGV == 3);


# should open files first, else
# if last file cannot open, it is
# a waste of time to read files before

open (IDPOS,  "<$ARGV[0]")
  or die "Cannot open ID_probe file: $!";
  
open (PRE, "<$ARGV[1]")
  or die "Cannot open preprocessed_File: $!";
  
open(OUT,">$ARGV[2]")
  or die "Cannot open output_File: $!";
 
 
print OUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup\n";
  
  my %chrs;
  my %starts;
  my %ends;
  my $line;
  my @parts;
  
  while($line = <IDPOS>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	$chrs{$parts[0]} = $parts[2];
	$starts{$parts[0]} = $parts[3];
	$ends{$parts[0]} = $parts[4];
  }

  close (IDPOS);
  
  while ($line = <PRE>)
  {
  	chomp $line;
	@parts = split /\t/, $line;
	if (exists $chrs{$parts[0]} )
	{
		print OUT "$chrs{$parts[0]}\t\.\t$parts[0]\t$starts{$parts[0]}\t$ends{$parts[0]}\t$parts[1]\t\.\t.\t.\n";			
	}
	else
	{
		next;
	}
  }
  close (PRE);
  close (OUT);
  print "\a";
  exit;
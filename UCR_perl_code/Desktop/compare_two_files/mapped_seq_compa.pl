#!/usr/bin/perl -w

# compare whether the  reads mapped
# to the cDNA and whole gene(include intron)

use strict;

my $usage = "$0 <cdna> <seq>";

die $usage unless(@ARGV == 2);

open (CDNA,  "<$ARGV[0]")
  or die "Cannot open former_File: $!";

open (SEQ, "<$ARGV[1]")
  or die "Cannot open redo_File: $!";
  
  my %gene;
  my %position;
  my $line;
  my @parts;
  my $same = 0;
  my $same_read_dif_gene = 0;
  my $same_read_dif_pos = 0;
  my $diff = 0;
  my $total_reads = 0;
  
  while($line = <CDNA>)
  {
	chomp $line;
	@parts = split /\t/, $line;
	$gene{$parts[0]} = $parts[7];
	$position{$parts[0]} = $parts[8];	
  }

  close (CDNA);
  
  while ($line = <SEQ>)
  {
  	$total_reads ++;
  	chomp $line;
	@parts = split /\t/, $line;
	if (exists $gene{$parts[0]} )
	{
		if ($gene{$parts[0] } eq $parts[7] && 
		     $position{$parts[0]}  == $parts[8])
		{
			$same ++;
		}
		else
		{
			if($gene{$parts[0] } ne $parts[7])
			{$same_read_dif_gene++; }
			else
			{ $same_read_dif_pos++}
		}	
	}
	else
	{
		$diff++;	
	}
  }
  
 
  close (SEQ);
  
   my $sum = $same+ $same_read_dif_gene +  $same_read_dif_pos +  $diff ;
  print "same map = $same, new_map = $diff \n ",
         "old map but to different gene = $same_read_dif_gene \t",
         "same but diff pos: $same_read_dif_pos\n",
         "sum = $sum\t total reads = $total_reads\n";
  print "\a";
  exit;
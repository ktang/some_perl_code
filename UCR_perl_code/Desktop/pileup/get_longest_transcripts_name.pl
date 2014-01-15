#!/usr/bin/perl -w

#input the transcripts-length file
# to filter the longest transcripts
# for each gene.

use strict;

my $usage = "$0 <transcript_length_file>";

die $usage unless(@ARGV == 1);

open (IN, '<', $ARGV[0])
  or die "Cannot open inFile: $!";

  my $line;
  my @parts;
  my %genes_length;
  my $gene;
  my %transcripts;
  
  while($line = <IN>)
  {
	chomp($line);
	 
	@parts = split /\t/,$line;

	$gene = substr($parts[0],0,9);
     	 
     	 if (exists $genes_length{$gene} )
     	 {
     	    if ($parts[1] > $genes_length{$gene})
     	    	{
     	    			$genes_length{$gene} = $parts[1];
     	    			$transcripts{$gene} = $parts[0];
     	    	}
     	 }
     	 else
     	 {
     	 	$genes_length{$gene} = $parts[1];
     	 	$transcripts{$gene} = $parts[0];;
     	 }
   }
  
  close (IN);
  
  foreach $gene(sort keys %transcripts)
  {
  	print $transcripts{$gene},"\n";
  }
  
  print "\a";
  exit;
  
  
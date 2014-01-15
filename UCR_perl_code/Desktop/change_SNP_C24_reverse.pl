#!/usr/bin/perl -w
# change old_C24_SNPs_in_cDNAs.txt the reverse stand
# to the + strand

use strict;

my $usage = "$0 <strand> <old> <new>";

die $usage unless(@ARGV == 3);

open (STRAND, "<$ARGV[0]")
  or die "Cannot open $ARGV[0]: $!";

open (IN, "<$ARGV[1]")
	or die "Cannot open $ARGV[1]:$!";

open (OUT,  ">$ARGV[2]")
  or die "Cannot open $ARGV[2]: $!";
  
  my @parts;
  my $line;
  my %strand;
  
  while($line = <STRAND>)
  {
	 chomp $line;
	 @parts = split /\t/,$line; 
    $strand{$parts[0]} = $parts[1];
  }
  
  close(STRAND);
  
  $line = <IN>;
  print OUT $line;
  
  while ($line = <IN>)
  {
  	 chomp $line;
	 @parts = split /\t/,$line;
	 
	 if (exists $strand{$parts[0]})
	 {
		 
	 	if($strand{$parts[0]} eq "+")
	 	{print OUT $line,"\n";}
	 
		 elsif ($strand{$parts[0]} eq "-")
  		{
			if($parts[2] eq "A")
			{$parts[2] = "T";} 
			
			elsif($parts[2] eq "T")
			{$parts[2] = "A";} 
			
			elsif($parts[2] eq "C")
			{$parts[2] = "G";} 
			
			elsif($parts[2] eq "G")
			{$parts[2] = "C";} 	
			
			else
			{print $line,"\n";}	
			
			if($parts[3] eq "A")
			{$parts[3] = "T";} 
			
			elsif($parts[3] eq "T")
			{$parts[3] = "A";} 
			
			elsif($parts[3] eq "C")
			{$parts[3] = "G";} 
			
			elsif($parts[3] eq "G")
			{$parts[3] = "C";} 	
			
			else
			{print $line,"\n";}	
			
			print OUT "$parts[0]\t$parts[1]\t$parts[2]\t$parts[3]\t$parts[4]\n";			
  		}
  		
  		else
  		{print $line,"\n";}	
  		
	}
	else	 {print $line,"\n";}	 	
 }

  close (IN);
  close (OUT);
  print "\a";
  exit;
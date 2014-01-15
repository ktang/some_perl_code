#!/usr/bin/perl -w
# input id_length for Col-0 and 
# C24, blastn result is the third 
# input, output should be 
# Col-0 id, length, and status
# 0: no hit
# 1: same
# 2: query name and subject name
#    are same, but seq are not the same
# 3: qurey name and subject name
#  are different
# NOTICE: for the length < 100, this 
# may be because pre-tRNA are the same.

use strict;


my $usage = "$0 <Col_0_length> <C24_length> <blastn_output>";
die $usage unless(@ARGV == 3);

open (IN1, '<', $ARGV[0])
  or die "Cannot open Col0_length: $!";
   
open (IN2, '<', $ARGV[1])
  or die "Cannot open C24_length: $!";

  my %C24_len;
  my %Col0_len;
  my %Col0_status;
  my $line;
  my @parts;
  
  
  while($line = <IN1>)
  {
  	chomp($line);
  	@parts = split /\t/,$line;
	$Col0_len{$parts[0]} = $parts[1] ;
	$Col0_status{$parts[0] } = 0;
  }
  
  close(IN1);
  
  
  
  while($line = <IN2>)
  {
     chomp($line);
     @parts = split /\t/,$line;
   	 $C24_len{$parts[0]} = $parts[1];
  }
	close (IN2);
 
 open (IN3,"<$ARGV[2]")
   or die "Cannot open blast_output:$!";
   
   while($line = <IN3>)
   {
   		chomp $line;
   		@parts = split /\t/, $line;
   		
   		if ( $parts[0] eq $parts[1])
   		{
   			if ($Col0_len{$parts[0]} == $C24_len{$parts[1]} 
   				&& $parts[2]  == 100 && $parts[6] == 1
   				&& $parts[7] == $Col0_len{$parts[0]}
   			    )	
   			    
   			{
   				$Col0_status{$parts[0]} = 1;
   			}
   			
   			elsif($Col0_status{$parts[0]} != 1 )
   			{
   					$Col0_status{$parts[0]} = 2;
   			}
   		}
   		
   		elsif ( $Col0_status{$parts[0]} == 0)
   		{
   			$Col0_status{$parts[0]} = 3;	
   		}
   	
   }
   
   close (IN3);
   
   foreach my $gene(sort keys %Col0_status)
   {
   		print "$gene\t$Col0_len{$gene}\t$Col0_status{$gene}\n"
   	}
 
  
  exit;
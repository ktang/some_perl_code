#!/usr/bin/perl -w

#using the exists transcripts and their
#position to extract form the original
#pileup file again

use strict;

my $usage = "$0 <gene_and_position> <pileup_in> ";
#<pileup_extracted>";
die $usage unless(@ARGV == 2);
my ($inFile1, $inFile2) = @ARGV[0..2];



open (IN1, "<$inFile1")
  or die "Cannot open $inFile1: $!";
  
open (IN2, "<$inFile2")
  or die "Cannot open $inFile2: $!";

open (OUT, ">$inFile2.May_8_extracted")
  or die "Cannot open outFile: $!";
  
  my @parts;
  my %genes;
  my $char;
  my $i;
  my $len;
  my $ratio;
  my $line;
  
 while($line = <IN1>)
 {
 	chomp($line);
 	@parts = split /\t/,$line;
 	$genes{$parts[0]} ->{$parts[1]} = -1;
 }
 
 while ($line = <IN2>)
 
 {
	chomp($line);
	
	
	@parts = split /\t/,$line;
	
	if (exists $genes{$parts[0]})
	{
		if (exists $genes{$parts[0]}-> {$parts[1]})
		{
	 		$len = length($parts[4]);
	 		my $count_match = 0;
			my ($num_A, $num_T, $num_C, $num_G) =(0,0,0,0);
			
			for ($i = 0 ; $i < $len; $i++)
    	  	{
     	 		$char = substr($parts[4], $i,1 );
   				if ( $char eq "\." ||  $char eq ",")
     	 		{
      				$count_match++;
      			}
      			elsif($char eq "a" || $char eq "A")
      			{
      				$num_A++;
      			}
      			elsif ($char eq "t" || $char eq "T")
      			{
      				$num_T++;
      			}
      		elsif ($char eq "c" || $char eq "C")
      		{
      			$num_C++;
      		}
      		elsif ($char eq "g" || $char eq "G")
      		{
      			$num_G++;
      		}
      	}
		
		
			 my $total = $count_match+$num_A+$num_T+$num_C+$num_G;
			 if ($total != $parts[3])
			 {
			 	print $line,"\n";
			 }
			 
			 if ($parts[3] >= 10 && ($count_match/$parts[3]) < 0.9 )
			 {
				print OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
				"\t",$parts[3],"\t",$parts[4],"\t",$count_match,
				"\t",$num_A,"\t",$num_T,"\t",
				$num_C,"\t",$num_G,"\t",$total,"\n";
			 }
		}
	 }
	 
}

print "\a";
exit;
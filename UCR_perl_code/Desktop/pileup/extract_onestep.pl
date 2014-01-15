#!/usr/bin/perl -w
#extract from the original pileup
#file, using the criteria that 
#(1) depth >= 20 $parts[3] >= 20
#(2) match percentage < 0.95

use strict;

my $usage = "$0 <original pileup>";
die $usage unless(@ARGV == 1);

open (IN, "<$ARGV[0]")
  or die "Cannot open inFile: $!";

open (OUT, ">$ARGV[0].onestep_extracted")
  or die "Cannot open outFile: $!";
  
  my @parts;
  my $len;
  my $char;
  my $i;
  my $count;
  my $ratio;
  my $line;
  
 while($line = <IN>)
 {
	chomp($line);
	$count = 0;
	
	@parts = split /\t/,$line;
	
	if ($parts[3] >= 20)
	{
	 	$len = length($parts[4]);
	 	$count = 0;
		my ($num_A, $num_T, $num_C, $num_G) =(0,0,0,0);
		for ($i = 0 ; $i < $len; $i++)
      	{
      		$char = substr($parts[4], $i,1 );
   			if ( $char eq "\." ||  $char eq ",")
     	 	{
      			$count++;
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
		
		$ratio = $count/$parts[3];
		
		if ( $ratio < 0.95)
		{
			 my $total = $count+$num_A+$num_T+$num_C+$num_G;
			 if ($total != $parts[3])
			 {
			 	print $line,"\n";
			 }
			 
			print OUT $parts[0],"\t",$parts[1],"\t",$parts[2],
		"\t",$parts[3],"\t",$parts[4],"\t",$count,
		"\t",$num_A,"\t",$num_T,"\t",
		$num_C,"\t",$num_G,"\t",$total,"\n";
		}
	 }
	 
}

print "\a";
exit;
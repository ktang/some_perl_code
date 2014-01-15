#!/usr/bin/perl -w

# caculate the freq for cut_table_rm0

use strict;

my $usage = "$0 <in> <out>";

die $usage unless(@ARGV == 2);


open (IN, "<$ARGV[0]")
  or die "Cannot open inFile: $!";

open (OUT, ">$ARGV[1]")
  or die "Cannot open outFile: $!";
  
  
  my $line;
  my @parts;
  
  $line = <IN>;
  
  print OUT "mRNAID	SNPPos	Col0Base	C24Base	WT(Col0XC24)_from_Col0	WT(Col0XC24)_from_C24	freq_WT(Col0XC24)_from_C24	WT(C24XCol0)_from_Col0	WT(C24XCol0)_from_C24	freq_WT(C24XCol0)_from_C24	nrpd1a(Col0XC24)_from_Col0	nrpd1a(Col0XC24)_from_C24	freq_nrpd1a(Col0XC24)_from_C24	nrpd1a(C24XCol0)_from_Col0	nrpd1a(C24XCol0)_from_C24	freq_nrpd1a(C24XCol0)_from_C24	nrpd1b(Col0XC24)_from_Col0	nrpd1b(Col0XC24)_from_C24	freq_nrpd1b(Col0XC24)_from_C24	nrpd1b(C24XCol0)_from_Col0	nrpd1b(C24XCol0)_from_C24	freq_nrpd1b(C24XCol0)_from_C24	drd1(Col0XC24)_from_Col0	drd1(Col0XC24)_from_C24	freq_drd1(Col0XC24)_from_C24	drd1(C24XCol0)_from_Col0	drd1(C24XCol0)_from_C24	freq_drd1(C24XCol0)_from_C24	total_depth\n";
  
  while($line = <IN>)
  {
	chomp($line);
	@parts = split /\t/,$line;
	print OUT "$parts[0]\t$parts[1]\t$parts[2]\t$parts[3]\t";
	
	for(my $i = 4; $i <=18; $i+=2)
	{
		if ($parts[$i]+ $parts[($i+1)]  == 0)
		{print OUT "0\t0\tNA\t";}
		
		else
		{
			my $ratio = $parts[($i+1)] /(  $parts[$i]+ $parts[($i+1)] );
			print OUT "$parts[$i]\t$parts[($i+1)]\t$ratio\t";
		}	
		
	}
 print OUT "$parts[20]\n";

 
}


  close (IN);
  close (OUT);
 
  exit;
#!/usr/bin/perl -w

use strict;


my $usage = "$0 <in1> <out>";
die $usage unless(@ARGV >= 2);
my ($in, $outFile) = @ARGV[0..1];

my $IN;

my $OUT;


open ($IN, '<', $in)
	or die "Cannot open $in:$!";


open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  


my $line;
my %genes;
my @parts;


while ($line = <$IN>)
{
	chomp($line);
	@parts = split /\t/,$line;
	$genes{$parts[0]}->{$parts[1]} = $parts[5]/$parts[3];
}


	foreach my $gene(sort keys %genes)
  	{
  		 my $flag = 0;
  		 my $flag_true = 1;
 		 my @k = keys %{$genes{$gene}};
 
 		if ($genes{$gene}->{$k[0]} < 0.2)
    	{
 		 	$flag = 0;
 		}
 		
 		elsif ($genes{$gene}->{$k[0]} < 0.4)
 	    {
 		 	$flag = 2;
 		 }
 		
 		elsif ($genes{$gene}->{$k[0]} < 0.6)
 		 {
 		 	$flag = 4;
 		 }
 		
 		elsif ($genes{$gene}->{$k[0]} < 0.8)
 		 {
 		 	$flag = 6;
 		 }
 		 
 		 elsif ($genes{$gene}->{$k[0]} <= 1)
 		 {
 		 	$flag = 8;
 		 }
 		 
 		 if ($flag == 0 ||$flag ==2 || $flag ==4 || $flag == 6 )
 		 {
 			 for (my $i = 1; $i <= $#k; $i++ )
 			 {
 		 		if (!( $genes{$gene}->{$k[$i]} >= ($flag*0.1) &&
 		 	       $genes{$gene}->{$k[$i]} < ( ($flag + 2)*0.1)))
 		    	{
 		    		$flag_true = 0;
 		    	}
 		 	}
 		 }
 		 
 		 elsif($flag == 8)
 		 {
 		 	for (my $i = 1; $i <= $#k; $i++ )
 			 {
 		 		if (!( $genes{$gene}->{$k[$i]} >= 0.8 &&
 		 	       $genes{$gene}->{$k[$i]} < 1 ))
 		    	{
 		    		$flag_true = 0;
 		    	}
 		 	}
 		 }
 		 
 		 if ($flag_true == 1)
 		 {
 		 	print $OUT $gene,"\t",$flag,"\n";
 		 }
 		 
 		 else
 		 {
 		 	print $gene,"\n";
 		 }
 	}


 
 print "\a";
 exit;
#!/usr/bin/perl -w

# count TE number in annotated file
#/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/result/send/label_ZDP_peaks_with_annotation.xls


use strict;

=head2
my $usage = "$0 <input1> <input2> <output>";

die $usage unless(@ARGV == 3);

open (IN1 , "<$ARGV[0]")
  or die "Cannot open input_file1: $!";
  
open (IN2,  "<$ARGV[1]")
  or die "Cannot open input_File2: $!";
  
open (OUT,">$ARGV[2]")
  or die "Cannot open output files: $!";
=cut  

open (IN , "/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/result/send/label_ZDP_peaks_with_annotation.xls")
  or die "Cannot open input_file: $!";
  
  my $l1;
 
  my @p1;
 
    
#  print OUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup\n";
 # $l1 = <IN1> ,$l2 = <IN2>;
  my ($total ,$te,$none) = (0,0,0); 
  my ($gene,$intergenic) = (0,0);
  
  <IN>;#chomp head
  while($l1 = <IN> )
  {
	chomp $l1;
	
	@p1 = split /\t/, $l1;
	$total++;
	
	if ($p1[7] eq "NONE") #TE = NONE, still can be TE
	
	{
		if($p1[5] eq "transposable_element_gene")
		{
			$te++;
		}
		
		
		else
		{
			$none ++;
		
			if ($p1[4] ne "NONE")
			{
				$gene++;
			#		$intergenic++;	
				print "gene:$l1\n";
			}	
		
			else
			{
					if(($p1[5] eq "NONE") and ($p1[9] ne "NONE") )
					{	
						$intergenic++;
					}
					else{print STDERR "error:$l1\n"; }
			}
		}
	}
	
	elsif ( $p1[7] =~ /AT.TE/ )
	{ 
		$te++;
	}
	
	else {print STDERR "error:$l1\n";}
}
  
  print STDERR "total = $total\nTE=$te\nNone=$none\ngene=$gene\ninter=$intergenic\n";

  close (IN);
  
    

  exit;
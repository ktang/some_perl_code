#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2> <in3> <out>";
die $usage unless(@ARGV >= 4);
my ($inFile1, $inFile2, $inFile3, $outFile) = @ARGV[0..3];

my $IN1;
my $IN2;
my $IN3;
my $OUT;

open ($IN1, '<', $inFile1)
  or die "Cannot open $inFile1: $!";
   
open ($IN2, '<', $inFile2)
  or die "Cannot open $inFile2: $!";
  
open ($IN3, '<', $inFile3)
  or die "Cannot open $inFile3: $!";

open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  
  my %genes;
  my %uniq;
  my %gene_depth;
  my $line;
  my @parts;
  my $i;
  my @name_len;
  
  while($line = <$IN1>)
  {
  	chomp($line);
  	@name_len = split /\t/,$line;
  	$genes{$name_len[0]} = $name_len[1];
  }
  
  while($line = <$IN2>)
  {
     
    @parts = split /\t/,$line;
    if (exists $genes{ $parts[8]})
    {
    		for($i = $parts[9]; $i < $parts[9] + $parts[6]; $i++ )
    		{
    			$gene_depth{ $parts[8]} ->{"first"} -> [$i] += $parts[1];
    		}
    	}
    }	
  
   while($line = <$IN3>)
  {
     
    @parts = split /\t/,$line;
    if (exists $genes{ $parts[8]})
    {
    		for($i = $parts[9]; $i < $parts[9] + $parts[6]; $i++ )
    		{
    			$gene_depth{ $parts[8]} -> {"second"} -> [$i] += $parts[1];
    		}
    	}
    }
  
  
  foreach my $gene(sort keys %genes)
  {
 	 print $OUT $gene,"\t",$genes{$gene},"\t";
 	 
		
		for($i = 0; $i < $genes{$gene}; $i++)
		{
		  if(!defined($gene_depth{$gene}->{"first"}->[$i]))
		  {
			$gene_depth{$gene} -> {"first"} -> [$i] = 0;
	   	  }
		  print $OUT $gene_depth{$gene}->{"first"}->[$i],"\t";
		 
		}
		 


		for($i = 0; $i < $genes{$gene} -1; $i++)
		{
		  if(!defined($gene_depth{$gene}->{"second"}->[$i]))
		  {
			$gene_depth{$gene}->{"second"}->[$i] = 0;
	   	  }
		  print $OUT $gene_depth{$gene}->{"second"}->[$i],"\t";
		 
		}
		 if(!defined($gene_depth{$gene}->{"second"}->[$genes{$gene} -1]))
		  {
			$gene_depth{$gene}->{"second"}->[$genes{$gene} -1] = 0;
	   	  }
		  print $OUT $gene_depth{$gene}->{"second"}->[$genes{$gene} -1],"\n";	   
 }
  
  
  
=pod 
  foreach my $gene(sort keys %genes)
  {
 	 print $OUT $gene,"\n";
 	 
		#my @size = $gene_depth{$gene};
		for($i = 0; $i < $genes{$gene}; $i++)
		{
		  print $OUT $i+1,"\t";
		}
		print $OUT "\n";
		for($i = 0; $i < $genes{$gene}; $i++)
		{
		if(!defined($gene_depth{$gene}->[$i]))
		{
		$gene_depth{$gene}->[$i] = 0;
		}
		  print $OUT $gene_depth{$gene}->[$i],"\t";
		}
	    print $OUT "\n";
 }
=cut

  close ($IN1);
  close ($IN2);
  close ($IN3);
  close ($OUT);
  print "\a";
  exit;
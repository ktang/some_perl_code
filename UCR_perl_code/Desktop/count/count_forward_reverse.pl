#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2> <out>";
die $usage unless(@ARGV >= 3);
my ($inFile1, $inFile2, $outFile) = @ARGV[0..2];

my $IN1;
my $IN2;
my $OUT;

open ($IN1, '<', $inFile1)
  or die "Cannot open $inFile1: $!";
   
open ($IN2, '<', $inFile2)
  or die "Cannot open $inFile2: $!";

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
 	   if (exists $uniq{$parts[8].$parts[0] } )
    	{
     	 next;
    	}
    	else
    	{
    		$uniq{$parts[8].$parts[0]} = -1;
    		if ($parts[7] eq "+")
    		{
    			$gene_map{ $parts[8]} -> {"+"} += $parts[1];
    		}
    		elsif ($parts[7] eq "-")
    		{
    			$gene_map{ $parts[8]} -> {"-"} += $parts[1];
    		}
    	}
    }	
  }
  
  
  foreach my $gene(sort keys %genes)
  { 
     if (!defined($gene_map{$gene} -> {"+"}))
     {
       $gene_map{$gene} -> {"+"} =0;
     }
     
     if (!defined($gene_map{$gene} -> {"-"}))
     {
       $gene_map{$gene} -> {"-"} =0;
     }
     
 	 print $OUT $gene,"\t",$gene_map{$gene} -> {"+"},"\t",$gene_map{$gene} -> {"-"},"\n";
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
  close ($OUT);
  print "\a";
  exit;
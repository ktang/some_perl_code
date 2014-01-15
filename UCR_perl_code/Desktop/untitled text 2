#!/usr/bin/perl -w
my $usage = "$0 <in1> <in2> <out>";
die $usage unless(@ARGV >= 2);
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
  my $gene_name;
  
  while($line = <$IN1>)
  {
  	chomp($line);
  	$genes{$line} = -1;
  }
  
  while($line = <$IN2>)
  {
     chomp($line);
     
    if ($line =~ m/^>/)
    { 
    	if ($line =~ m/(AT[1-5|M|C]G\d+)\.\d/)
   	 	{
		   	$gene_name = $1;		      
		}
   		
   		if (exists $genes{$gene_name})
   		{
   		 	if($line =~ m/chr[1-5|M|C]:(\d+)-(\d+)/ )
   		 	{
	   		 	$genes{$gene_name} = $2 - $1 + 1;
	   		}   		 	
   		}    	
    }  
  }


 foreach my $gene(sort keys %genes)
  {
 	 
	    print $OUT $gene,"\t",$genes{$gene},"\n";
	  
  }

  close ($IN1);
  close ($IN2);
  close ($OUT);
  print "\a";
  exit;
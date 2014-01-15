#!/usr/bin/perl -w
#with the gene name in <gene_name>
#to get all the position for all the
#1236 genes, inDir contains all the extracted 
#position from the original pileup file
#(depth >= 20 and match percentage < 95%

use strict;

my $usage = "$0 <gene_name> <inDir> <outFile>";
die $usage unless(@ARGV >= 3);
my ($infile, $indir, $outfile) = @ARGV[0..2];

my $IN;
my $OUT;

open($IN, '<' ,$infile)
	 or die "Cannot open $infile: $!";

open ($OUT, '>', $outfile)
  or die "Cannot open $outfile: $!";
  
  my @parts;
  my %genes;
  my $line;
  my $count = 0;
  
  
 while($line = <$IN>)
 {
  	chomp($line);
  	$genes{$line} = -1;
 }
 
 close($IN);
 
 opendir (INDIR, $indir)
 	or die "Cannot open dir $indir:$!";
 	
 my @files = grep {/extracted$/} readdir INDIR;
 	
 foreach my $file(@files)
 {
 	open ($IN, '<', "$indir/$file")
 		or die "Cannot open $file: $!";
 	
 	while ($line = <$IN>)
 	{
 		chomp($line);
		@parts = split /\t/,$line;
		my $name = substr($parts[0], 0, 9);
		
		
		if (exists $genes{$name} )
		{
			if($genes{$name} == -1)
			{
				$count++;
				$genes{$name} = 1;
			}
		
			print $OUT $parts[0],"\t",$parts[1],"\n";
			
		}
	}
	close($IN);
}
 
 print "total detected genes:$count\n";
 close ($OUT);

print "\a";
exit;
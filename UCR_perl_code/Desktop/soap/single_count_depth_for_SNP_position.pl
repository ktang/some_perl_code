#!/usr/bin/perl -w
# count the depth of SNP position 
# directly from the unique match
# soap output

use strict;

my $usage = "$0 <soap_output> <SNP_depth_output>";

die $usage unless(@ARGV == 2);

my %genes;
my $line;
my @parts;

open (IN,"<$ARGV[0]") or die "Cannot open the input soapfile $ARGV[0]:$!";

open (OUT, ">$ARGV[1]") or die "Cannot opne the output file $ARGV[1]:$!";

open(SNP,"</mnt/disk4/genomes/heterosis/kai/info/new_C24_SNPs_in_cDNAs.txt") or die "cannot open SNP file:$!";

<SNP>;

while ($line = <SNP>)
{
	chomp $line;
	@parts = split /\t/,$line;
	$genes{$parts[0]} -> {$parts[1]} = 0;
}

close (SNP);

while ($line = <IN>)
{
	chomp $line;
	@parts = split /\t/, $line;
	my $start = $parts[8];
	my $end = $start + $parts[5] -1;
	if (exists $genes{$parts[7]})
	{
				
		
		foreach my $position(sort keys %{$genes{$parts[7]}})
		{
			
			if ($position >= $start && $position <= $end )
			{
				($genes{$parts[7]}->{$position})++;
				
			}
		}	
	}
}


close(IN);

foreach my $gene(sort keys %genes)
{
	foreach my $position(sort {$a<=>$b} keys %{$genes{$gene}})
	{
		print OUT "$gene\t$position\t$genes{$gene}->{$position}\n";	
	}
}

close (OUT);
print "\a";
exit;
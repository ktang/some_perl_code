#!/usr/bin/perl -w
# count the depth of SNP position 
# directly from the unique match
# soap output

use strict;

my $usage = "$0 <SNP_input> <output_sorted>";

die $usage unless(@ARGV == 2);

my %genes;
my $line;
my @parts;

open (IN,"<$ARGV[0]") or die "Cannot open the input soapfile $ARGV[0]:$!";

open (OUT, ">$ARGV[1]") or die "Cannot opne the output file $ARGV[1]:$!";



<IN>;

while ($line = <IN>)
{
	chomp $line;
	@parts = split /\t/,$line;
	$genes{$parts[0]} -> {$parts[1]} = "$parts[2]\t$parts[3]\t$parts[4]";
}

close (IN);

print OUT "mRNAID\tSNPPos\tCol0Base\tC24Base\tSNPID\n";

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
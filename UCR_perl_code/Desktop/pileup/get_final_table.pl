#!/usr/bin/perl -w

# use C24_SNPs_in_cDNAs.txt and pos_lib_info_May17.txt
# to get the final result 

use strict;


my $usage = "$0 <C24_SNPs_in_cDNAs.txt> <pos_lib_info_May17.txt> <output>";

die $usage unless(@ARGV == 3);

my $line;

my %c24;
my %col0;
my @parts;

open (INSNP, "<$ARGV[0]")
	or die "Cannot open $ARGV[0]:$!";

while ($line = <INSNP>)
{
	chomp($line);
	@parts = split /\t/,$line;
	$col0{$parts[0].$parts[1]} = $parts[2];
	$c24{$parts[0].$parts[1]} = $parts[3];
}

close(INSNP);

open (IN,"<$ARGV[1]") or die "cannot open $ARGV[1]:$!";

open (OUT, ">$ARGV[2]")
  or die "Cannot open $ARGV[2]: $!";
  
print OUT "gene\tposition\tCol0_base\tC24_base\tdepth_Col0_flowcell54_lane1\tdepth_C24_flowcell54_lane1\tfreq_C24_flowcell54_lane1\tdepth_Col0_flowcell54_lane2\tdepth_C24_flowcell54_lane2\tfreq_C24_flowcell54_lane2\tdepth_Col0_flowcell54_lane3\tdepth_C24_flowcell54_lane3\tfreq_C24_flowcell54_lane3\tdepth_Col0_flowcell54_lane5\tdepth_C24_flowcell54_lane5\tfreq_C24_flowcell54_lane5\tdepth_Col0_flowcell54_lane6\tdepth_C24_flowcell54_lane6\tfreq_C24_flowcell54_lane6\tdepth_Col0_flowcell54_lane7\tdepth_C24_flowcell54_lane7\tfreq_C24_flowcell54_lane7\tdepth_Col0_flowcell54_lane8\tdepth_C24_flowcell54_lane8\tfreq_C24_flowcell54_lane8\tdepth_Col0_flowcell55_lane1\tdepth_C24_flowcell55_lane1\tfreq_C24_flowcell55_lane1\n";
  
while($line = <IN>)
{
	chomp $line;
	@parts = split /\t/,$line;
	if (exists $col0{$parts[0].$parts[1]} )
	{
		print OUT "$parts[0]\t$parts[1]\t$col0{$parts[0].$parts[1]}\t$c24{$parts[0].$parts[1]}\t";
		
		if ($parts[2]+$parts[3] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[2]\t$parts[3]\t",$parts[3]/($parts[2]+$parts[3]),"\t";}
		
		if ($parts[4]+$parts[5] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[4]\t$parts[5]\t",$parts[5]/($parts[4]+$parts[5]),"\t";}
		
		if ($parts[6]+$parts[7] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[6]\t$parts[7]\t",$parts[7]/($parts[6]+$parts[7]),"\t";}
		
		if ($parts[8]+$parts[9] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[8]\t$parts[9]\t",$parts[9]/($parts[8]+$parts[9]),"\t";}
		
		if ($parts[10]+$parts[11] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[10]\t$parts[11]\t",$parts[11]/($parts[10]+$parts[11]),"\t";}
		
		if ($parts[12]+$parts[13] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[12]\t$parts[13]\t",$parts[13]/($parts[12]+$parts[13]),"\t";}
		
		if ($parts[14]+$parts[15] == 0)
		{print OUT "0\t0\tNA\t";}
		else
		{print OUT "$parts[14]\t$parts[15]\t",$parts[15]/($parts[14]+$parts[15]),"\t";}
		
		if ($parts[16]+$parts[17] == 0)
		{print OUT "0\t0\tNA\n";}
		else
		{print OUT "$parts[16]\t$parts[17]\t",$parts[17]/($parts[16]+$parts[17]),"\n";}
		
		
		
	}
	
	else
	{print $line,"\n";}
	
}


close (IN);
close(OUT);
 
 print "\a";
 exit;
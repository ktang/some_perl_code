#!/usr/bin/perl -w
#use the C24_SNPs_in_cDNAs.txt to 
# get all the SNP position in the 
#original pileup file

use strict;

my $usage = "$0 <C24_SNP> <inDir> <outDir>";
die $usage unless(@ARGV == 3);

open (IN, "<$ARGV[0]")
  or die "Cannot open C24_SNP: $!";

my %positions;
my @parts;
my $line;
  
 while($line = <IN>)
 {
	chomp($line);
	@parts = split /\t/,$line;
	$positions{$parts[0].$parts[1]} = -1;
 }
 
 close(IN);
 
opendir (INDIR, $ARGV[1]) or die "Cannot open dir $ARGV[1]:$!";

my @files = grep {/\.pileup$/} readdir INDIR;
foreach my $file(@files)
{
     if($file =~ /(\S+)\.pileup$/)
     {
		 my $pre = $1;
		 my $output = $pre . ".pileup_SNP";
		 open (INFILE,"<$ARGV[1]/$file") or die "cannot open $file:$!";
		 open (OUTFILE, ">$ARGV[2]/$output") or die "cannot opne $output:$!";
		
		while($line = <INFILE>)
		{
			chomp $line;
			@parts = split /\t/,$line;
			if (exists $positions{$parts[0].$parts[1]}	)
			{
				print OUTFILE $line,"\n";	
			}
		}		 
		 
		 close(INFILE);
		 close(OUTFILE);
		 
	 }
}

print "\a";
exit;
#!/usr/bin/perl -w
# batch calculate each position 
# the match percentage

use strict;

my $usage = "$0 <inDir> <outDir>";
die $usage unless(@ARGV == 2);
my ($indir, $outdir) = @ARGV[0..1];

my $IN;
my $OUT;

my $line;
 
 opendir (INDIR, $indir)
 	or die "Cannot open dir $indir:$!";
 	
 my @files = grep {/^fl/} readdir INDIR;
 	
 foreach my $file(@files)
 {
 	my $outfile = $file."_May8_percentage";
 	open ($IN, '<', "$indir/$file")
 		or die "Cannot open $file: $!";
 	
 	open ($OUT, '>', "$outdir/$outfile")
 		or die "Cannot open $outfile: $!";
 	
 	while ($line = <$IN>)
 	{
 		chomp $line;
 		my @parts = split /\t/,$line;
 		my $ratio = $parts[5]/$parts[3];
 		print $OUT $line,"\t",$ratio,"\n";
 	
 	}
	
	close($IN);
	close ($OUT);
}
 

print "\a";
exit;
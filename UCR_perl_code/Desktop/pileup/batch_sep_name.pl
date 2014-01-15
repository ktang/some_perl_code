#!/usr/bin/perl -w

use strict;

my $usage = "$0 <inDir> <outDir>";
die $usage unless(@ARGV >= 2);
my ($indir, $outdir) = @ARGV[0..1];

my $IN;
my $OUT;

my $line;
 
 opendir (INDIR, $indir)
 	or die "Cannot open dir $indir:$!";
 	
 my @files = grep {/^class/} readdir INDIR;
 	
 foreach my $file(@files)
 {
 	my $outfile = $file."_sep";
 	open ($IN, '<', "$indir/$file")
 		or die "Cannot open $file: $!";
 	
 	open ($OUT, '>', "$outdir/$outfile")
 		or die "Cannot open $outfile: $!";
 	
 	while ($line = <$IN>)
 	{
 		$line =~ s/\./\t/;
 		 print $OUT $line;
	}
	close($IN);
	 close ($OUT);
}
 

print "\a";
exit;
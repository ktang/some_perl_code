#!/usr/bin/perl -w

# use the sorted probe file to get sorted gff file.

use strict;

my $usage = "$0 <in_dir> <out_dir>";
die $usage unless(@ARGV == 2);
my ($indir, $outdir) = @ARGV[0..1];

open (SORTED,"/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/sorted_position_ID_v2")
or die "cannot open sorted probe file:$!";

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.sgr$/} readdir INDIR;

my $line;
my %chrs;
my %starts;
my %ends;
my @pts;

my @sorted = <SORTED>;
close(SORTED);

foreach my $file(@files)
{
	#my %vals;
	my $i = 0;
	if($file =~ /(\S+)\.sgr$/)
	{
		 my $pre = $1;
		 my $output = $pre . ".gff";
		 
		 open(OUT, ">$outdir/$output")
		   or die "cannot open output file $output:$!";
		   
		 open (IN, "$indir/$file")
		   or die "cannot open input file $file:$!";
		   
		 while ($line = <IN>)
		 {
		 	chomp $line;
		 	@pts = split /\t/,$line;
		 	my $thisProbe = $sorted[$i];
		 	$i++;
		 	chomp $thisProbe;
		 	my @pts_probe = split "\t" , $thisProbe;
		 	if ( ($pts[0] eq $pts_probe[0] )and ($pts[1] eq $pts_probe[1]))
		 	{
		 		print OUT "$pts[0]\t.\t.\t$pts[1]\t$pts_probe[2]\t$pts[2]\t.\t.\t.\n";
		 	}	
		 	else
		 	{die "wrong order"}
		 }	
	
		close(IN);
		close(OUT);
	} 
}

print STDERR "\a";
exit;
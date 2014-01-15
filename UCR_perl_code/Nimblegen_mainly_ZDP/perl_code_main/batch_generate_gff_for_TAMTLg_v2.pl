#!/usr/bin/perl -w

# use the sorted probe file to get sorted gff file.

use strict;

my $usage = "$0 <in_dir> <out_dir>";
die $usage unless(@ARGV == 2);
my ($indir, $outdir) = @ARGV[0..1];

open (SORTED,"/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/sorted_position_ID_v2")
or die "cannot open sorted probe file:$!";

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/filter\.txt$/} readdir INDIR;

my $line;
my %chrs;
my %starts;
my %ends;
my @pts;

my @sorted = <SORTED>;
close(SORTED);

foreach my $file(@files)
{
	my %vals;
	if($file =~ /(\S+)\.txt$/)
	{
		 my $pre = $1;
		 my $output = $pre . ".gff";
		 
		 open(OUT, ">$ARGV[1]/$output")
		   or die "cannot open output file $output:$!";
		   
		 open (IN, "<$ARGV[0]/$file")
		   or die "cannot open input file $file:$!";
		   
		 while ($line = <IN>)
		 {
		 	chomp $line;
		 	@pts = split /\t/,$line;
		 	$vals{$pts[0]} = $pts[1];
		 }	
	
		close(IN);

		foreach my $thisline (@sorted)
		{
			chomp $thisline;
			@pts = split "\t", $thisline;
			if( exists $vals{$pts[3]})	
			{
				print OUT "$pts[0]\t.\t.\t$pts[1]\t$pts[2]\t$vals{$pts[3]}\t$pts[3]\t.\t.\n";
			}
			else 
			{
				print "ERROR: $thisline\n";
			}
		}	 	
	 
		close(OUT);
	} 
}

print STDERR "\a";
exit;

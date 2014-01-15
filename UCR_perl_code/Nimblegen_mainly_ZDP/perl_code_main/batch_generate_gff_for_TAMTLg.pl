#!/usr/bin/perl -w
# batch convert fastq file to fasta file
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
my @parts;


while ($line = <SORTED>)
{
	chomp $line;
	@parts = split /\t/,$line;
	$chrs{$parts[3]} = $parts[0];
	$starts{$parts[3]} = $parts[1];
	$ends{$parts[3]} = $parts[2];
 	}

foreach my $file(@files){
     if($file =~ /(\S+)\.txt$/){
		 my $pre = $1;
		 my $output = $pre . ".gff";
		 
		 open(OUT, ">$ARGV[1]/$output")
		   or die "cannot open output file $output:$!";
		   
		 open (IN, "<$ARGV[0]/$file")
		   or die "cannot open input file $file:$!";
		   
		 while ($line = <IN>)
		 {
		 	chomp $line;
		 	@parts = split /\t/,$line;
		 	if (exists $chrs{$parts[0]})
		 	{
		 		print OUT "$chrs{$parts[0]}\t.\t.\t$starts{$parts[0]}\t$ends{$parts[0]}\t$parts[1]\t.\t.\t$parts[0]\n";
		 	}
		 }
		 
	 
	 
	 close(IN);
	 close(OUT);
	 }
}

print "\a";
exit;

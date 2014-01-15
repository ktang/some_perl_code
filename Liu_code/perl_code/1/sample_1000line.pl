#!/usr/bin/perl -w

# get unique seqs
use strict;
use File::Spec;
use Bio::SeqIO;

my $usage = "$0 <input dir> <output dir>";
die $usage unless(@ARGV >= 2);
my ($inDir, $outDir) = @ARGV[0..1];

opendir(IND, $inDir) or die "Cannot open $inDir: $!";
my @inFiles = grep {/\.fasta/} readdir IND;
foreach my $inFile(@inFiles){

	my @parts = split /\./, $inFile;
	my $ext = pop @parts;
	my $outFile = "$outDir$parts[0]_first1000\.fasta";
	my $input = "$inDir/$inFile";
        my $OUT;
        my $IN;
        open($OUT, '>', $outFile)
            or die "cannot open $outFile:$!";
        open($IN, '<', $input)
            or die "cannot open $inFile:$!";
        for (my $i = 0; $i < 20; $i++)
        {
           my $line = <$IN>;
            print $OUT "$line";
        }
        close ($OUT);
    }

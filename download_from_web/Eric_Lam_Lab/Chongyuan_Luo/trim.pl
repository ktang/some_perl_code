#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $input = {};   

if (scalar(@ARGV) <= 0) {
    option();
}

else
{
my $m=""; my $d="";
GetOptions ('m=s'=>\$m, 'd=s'=>\$d);
my $dirlength=length $d; my $lastchar=substr($d, $dirlength-1, 1);
if ($lastchar cmp "/") {$d=$d."/";}

my $switch; my $genecount=0; my $max=0; my @read; my @gene; my $length; my %gene;
for (my $file=1; $file<=5; $file++)
{
    print "Loading Annotation of Chromosome $file...\n";
    open at_anno, "/Users/chongyuanluo/lab/Arabidopsis Genome/Tair 8/Annotation/NCBI_chr$file.tbl" or die $!;
    $switch=0;
    while (<at_anno>)
    {
        chop $_;
        @read=split("\t", $_);
        if (($read[2])and($read[2]eq"gene"))
        {
            $switch=1;
            $gene[$genecount][0]=$file;
            $gene[$genecount][1]=$read[0];$gene[$genecount][2]=$read[1];
        }
        elsif ($switch==1)
        {
            $switch=0; $gene[$genecount][3]=$read[4];
            if ($gene[$genecount][2]>$gene[$genecount][1])
            {
                $gene[$genecount][4]=$gene[$genecount][2]-$gene[$genecount][1];
            }
            else
            {
                $gene[$genecount][4]=$gene[$genecount][1]-$gene[$genecount][2];
            }
            $gene{$gene[$genecount][3]}=$gene[$genecount][4];
            $genecount++;
        }
    }
    close at_anno;
    print "$genecount\tDONE\n";
}
$genecount--;

my $outputname=$d.$m.".trim";
open file_out, ">$outputname" or die $!;
open file_in, "$m" or die $!;
$_=<file_in>; my $cutoff;
print file_out "$_";
while (<file_in>)
{
    chop $_;
    @read=split(/\t/,$_);
    if ($read[1]eq"")
    {
        print file_out "$_\n";    
    }
    else
    {
        $cutoff=int($gene{$read[0]}/50)+20;
        if ($cutoff<120)
        {
            for (my $run=$cutoff; $run<=120; $run++)
            {   $read[$run]="NA"; }
        }
        for (my $run=0; $run<=120; $run++)
        {
            print file_out "$read[$run]";
            print file_out "\t" if ($run<120);
        }
        print file_out "\n";
    }
}
close file_in;
close file_out;
}

sub option
{
    print "Trim.pl\nChongyuan Luo, Biotech, Rutgers\nFeb. 2010\n\nOptions\n";
    print "-m input pattern file\n";
    print "-d output path\n";
}
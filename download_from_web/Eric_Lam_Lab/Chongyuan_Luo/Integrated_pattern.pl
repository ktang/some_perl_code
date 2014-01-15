#!/usr/bin/perl -w
use strict;
my @read; my %mg;
open masked_gene, "/Users/chongyuanluo/lab/Arabidopsis Genome/RepeatMasked TIGR5/araTha5.tab.overlap" or die $!;
$_=<masked_gene>;
while (<masked_gene>)
{
    chop $_;
    @read=split(/\t/,$_);
    $mg{$read[0]}=1;
}
close masked_gene;

my $switch; my %gene; my $genecount; my @coverage;
for (my $file=1; $file<=5; $file++)
{
    print "Loading Annotation of Chromosome $file...\n";
    open at_anno, "/Users/chongyuanluo/lab/Arabidopsis Genome/Tair 8/Annotation/NCBI_chr$file.tbl" or die $!;
    $switch=0; $genecount=0;
        
        while (<at_anno>)
        {
            chop $_;
            @read=split("\t", $_);
            if (($read[2])and($read[2]eq"gene"))
            {
                $switch=1;
            }
            elsif ($switch==1)
            {
                $switch=0;
                for (my $run=0; $run<=8; $run++)
                { $gene{$read[4]}[$run]=0; }
                $genecount++;
            }
        }
        $genecount--;
}
print "DONE\n";

my @list;
@list=("H3K4me2_peaks.xls.overlap","H3K4me3_peaks.xls.overlap","H3K9me2_peaks.xls.overlap","H3K9Ac_peaks.xls.overlap",
       "H3K27me1_peaks.xls.overlap","H3K27me3_peaks.xls.overlap","H3K36me2_peaks.xls.overlap","H3K36me3_peaks.xls.overlap",
       "5mC_domain.txt.overlap");

for (my $run=0; $run<=8; $run++)
{
    open in_gene, "$list[$run]" or die $!;
    $_=<in_gene>;
    while (<in_gene>)
    {
        chop $_;
        @read=split(/\t/, $_);
        $gene{$read[0]}[$run]=1;
    }
    close in_gene;
}

open out_pattern, ">combined_pattern.txt" or die $!;
print out_pattern "\tH3K4me2\tH3K4me3\tH3K9me2\tH3K9Ac\tH3K27me1\tH3K27me3\tH3K36me2\tH3K36me3\t5mC\n";
foreach (keys %gene)
{
    if (! exists $mg{$_})
    {
        print out_pattern "$_\t";
        for (my $run=0; $run<=8; $run++)
        {
            print out_pattern "$gene{$_}[$run]";
            print out_pattern "\t" if ($run<8);
        }
        print out_pattern "\n";
    }
}
close out_pattern;
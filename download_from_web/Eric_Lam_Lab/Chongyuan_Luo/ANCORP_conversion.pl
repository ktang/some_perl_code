#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw[min max];

my $input = {};   

if (scalar(@ARGV) <= 0) {
    option();
}

else
{
my $r=""; my $c=""; my $n=""; my $a=""; my $d=""; my $l="";
GetOptions ('r=s'=>\$r, 'c=s'=>\$c, 'n=s'=>\$n, 'a=s'=>\$a, 'd=s'=>\$d, 'l=s'=>\$l);
my $dirlength=length $d; my $lastchar=substr($d, $dirlength-1, 1);
if ($lastchar cmp "/") {$d=$d."/";}

open ref_cdt, "$r" or die $!;
open input_cdt, "$c" or die $!;
my $outputname=$n."_".$a.".vis";
open out_cdt, ">$outputname" or die $!;

my @read; my $switch; my %gene; my $length;
if ($l==1)
{
    print "Loading Arabidopsis Genome Annotation...\n";
    for (my $file=1; $file<=5; $file++)
    {
        open at_anno, "/Users/chongyuanluo/lab/Arabidopsis Genome/Tair 8/Annotation/NCBI_chr$file.tbl" or die $!;
        $switch=0;
            while (<at_anno>)
            {
                chop $_;
                @read=split("\t", $_);
                if (($read[2])and($read[2]eq"gene"))
                {$switch=1; $length=max($read[0],$read[1])-min($read[0],$read[1]);}
                elsif ($switch==1) {$switch=0;$gene{$read[4]}=$length;}
            }
        close at_anno;
    }
    print "DONE\n";
}

$_=<input_cdt>;
chop $_;
$_="0".$_ if (substr($_,0,1) eq "\t");
print out_cdt "$_\n";

$_=<ref_cdt>;
print "Loading Reference pattern file...\n";
my @order; my $count=0;
while (<ref_cdt>)
{
    chop $_;
    @read=split(/\t/, $_);
    $order[$count]=$read[0];
    $count++;
}
$count--;
print "DONE\n";

print "Loading Input pattern file...\n";
my %input; 
while (<input_cdt>)
{
    chop $_;
    @read=split(/\t/, $_);
    $length=@read;
    $read[0]=~tr/gt/GT/;
    if ($l==1)
    {
        my $bin_length;
        $bin_length=int($gene{$read[0]}/50)+20;
        $bin_length=120 if ($bin_length>120);
        $read[$bin_length]=-50;
    }
    $input{$read[0]}=[@read];
}
print "DONE\n";

print "Performing Conversion...\n";
foreach (@order)
{
    if (exists $input{$_})
    {
        for (my $run=0; $run<=$length-1; $run++)
        {
            print out_cdt "$input{$_}->[$run]";
            print out_cdt "\t" if ($run<$length-1);
        }
        print out_cdt "\n";
    }
    else
    {
        print out_cdt "$_\t";
        for (my $run=1; $run<=$length-1; $run++)
        {
            print out_cdt "\t" if ($run<$length-1);
        }
        print out_cdt "\n";
    }

}
print "DONE\n";
}

sub option
{
    print "ANCORP_conversion.pl\nChongyuan Luo, Biotech, Rutgers\nSep. 2009\n\nOptions\n";
    print "-r Reference pattern file\n";
    print "-c Correlative pattern file\n";
    print "-n Name of the correlative pattern\n";
    print "-a Name of the anchor pattern\n";
    print "-l 1: Mark the TTS of transcription units\n";
    print "-d output path\n";
}
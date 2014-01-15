#!/usr/bin/perl -w
use strict;
use List::Util qw[min max];

use strict;
use Getopt::Long;

my $input = {};   

if (scalar(@ARGV) <= 0) {
    option();
}

else

{

my $m; my $d=""; my $b="";
GetOptions ('m=s'=>\$m, 'd=s'=>\$d, 'b=s'=>\$b);
my $dirlength=length $d; my $lastchar=substr($d, $dirlength-1, 1);
if ($lastchar cmp "/") {$d=$d."/";}

my $outputname=$d.$m.".overlap";

my @read; my $switch; my @gene; my $count=0; my @chromosome;

print "\nLoading Arabidopsis Genome Annotation...\n";
for (my $file=1; $file<=5; $file++)
{
    $chromosome[$file]=$count;
    open at_anno, "/Users/chongyuanluo/lab/Arabidopsis Genome/Tair 8/Annotation/NCBI_chr$file.tbl" or die $!;
    $switch=0;
        while (<at_anno>)
        {
            chop $_;
            @read=split("\t", $_);
            if (($read[2])and($read[2]eq"gene"))
            {$switch=1;$gene[$count][0]=$file;$gene[$count][1]=min($read[0],$read[1]);$gene[$count][2]=max($read[0],$read[1]);}
            elsif ($switch==1) {$switch=0; $gene[$count][3]=$read[4]; $gene[$count][4]=0; $count++;}
        }
    close at_anno;
}
print "DONE\n";
$chromosome[6]=$count;
$count--;

open peak_file, "$m" or die $!;

my $line=-1; my $chr; my $start; my $end; my $scan; my $found; my $quit; my $positive=-1; my $fc;
while (<peak_file>)
{
    $line++;
    if ($line>=18)
    {
        chop $_;
        @read=split("\t", $_);
        $chr=substr($read[0],3,1);
        $start=$read[1]; $end=$read[2]; $fc=$read[7];
        
        $scan=$chromosome[$chr]; $quit=0;
        while (($scan<$chromosome[$chr+1])and($quit==0))
        {
            if ($end<$gene[$scan][1]-$b) {$quit=1;}
            elsif (!($start>$gene[$scan][2]+$b))
            {
                $positive++ if ($gene[$scan][4]==0);
                $gene[$scan][4]=1;
            }
            $scan++;
        }
    }
}
close peak_file;

open gene_call, "> $outputname" or die $!;
print gene_call "$positive\n";
for (my $run=0; $run<=$count; $run++)
{
    if ($gene[$run][4]>0)
    {
        print gene_call "$gene[$run][3]\n";
    }
}

close gene_call;

}
sub option
{
    print "Peaks Alignment\nChongyuan Luo, Biotech, Rutgers\nFeb. 2009\n\nOptions\n";
    print "-m input peak files exported by MACS\n";
    print "-b maximum distance between peak and gene allowed\n";
    print "-d output path\n";
}
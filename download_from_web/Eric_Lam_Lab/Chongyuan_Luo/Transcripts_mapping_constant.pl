#!/usr/bin/perl -w
#Transcripts_mapping_constant.pl
#Sep. 2009
# Chongyuan Luo, Biotech, Rutgers
use strict;
use Getopt::Long;

my $input = {};   

if (scalar(@ARGV) <= 0) {
    option();
}

else
{
my $f=""; my $r=""; my $d="";
GetOptions ('f=s'=>\$f, 'r=s'=>\$r, 'd=s'=>\$d);
my $dirlength=length $d; my $lastchar=substr($d, $dirlength-1, 1);
if ($lastchar cmp "/") {$d=$d."/";}

my @read; my $switch; my @gene; my $genecount; my @coverage; my $min=9999;
my @head; 

my @out_pattern;
my @outputname=($d.$f.".constant",$d.$r.".constant");
open $out_pattern[1], ">$outputname[0]" or die $!;
open $out_pattern[0], ">$outputname[1]" or die $!;
foreach (@out_pattern)
{
    print $_ "\t";
    for (my $run=1; $run<=120; $run++)
    {
        print $_ "$run";
        print $_ "\t" if ($run<120);
    }
    print $_ "\n";
}

my $strand;
for (my $cycle=0; $cycle<=1; $cycle++)
{
    $strand=1 if ($cycle==0);
    $strand=-1 if ($cycle==1);
    
    open cov_chr, "$f" or die $! if ($strand==1);
    open cov_chr, "$r" or die $! if ($strand==-1);
    
    for (my $file=1; $file<=5; $file++)
    {
        print "Loading Annotation of Chromosome $file...\n";
        open at_anno, "/raid/home/chongyuan/Reference/NCBI_chr$file.tbl" or die $!;
        $switch=0; $genecount=0;
            
            while (<at_anno>)
            {
                chop $_;
                @read=split("\t", $_);
                if (($read[2])and($read[2]eq"gene"))
                {
                    $switch=1;
                    if ($read[0]<$read[1]) {$gene[$genecount][0]=1;} else {$gene[$genecount][0]=-1;}
                    $gene[$genecount][1]=$read[0];$gene[$genecount][2]=$read[1];
                }
                elsif ($switch==1)
                {
                    $switch=0; $gene[$genecount][3]=$read[4];
                    $genecount++;
                }
            }
            $genecount--;
            
        close at_anno;
        print "$genecount\tDONE\n"; 
        
        print "Loading coverage map for Chr$file...\n";
        my $basecount=-1; my $read=1;
        while ($read==1)
        {
            $_=<cov_chr>;
            $read=0 if (!$_);
            if ($read==1)
            {
                $read=0 if (substr($_,0,1) eq "S");
                if ($read==1)
                {
                    chop $_;
                    $basecount++;
                    @read=split(/\t/,$_);
                    $coverage[$basecount]=$read[1];
                }
            }        
        }
        print "$basecount\n";
        print "DONE\n";
        
        print "Calculating Coverage Over Gene Bodies $file...\n";
        my $head; my $tail; my $step; my $length; my $ratio; my $left; my $store; my $sum; my $baseleft; my $pos;
        
        for (my $gene=0; $gene<=$genecount; $gene++)
        {
            $head=$gene[$gene][1]-($gene[$gene][0]*1000); $tail=$gene[$gene][1]+($gene[$gene][0]*5000);
            if (($head<$basecount)and($head>0)and($tail>0)and($tail<$basecount))
            {
                for (my $run=0; $run<=120; $run++) {$head[$run]=0;}
            
                $step=1;
                while ($step<=5999)
                {
                    $pos=int($step/50);
                    $head[$pos]=$head[$pos]+($coverage[$head]/50);
                    $head=$head+$gene[$gene][0];
                    $step++;
                }             
                
                my $out_strand;
                $out_strand=$gene[$gene][0]*$strand;
                $out_strand=0 if ($out_strand==-1);
                
                print {$out_pattern[$out_strand]} "$gene[$gene][3]\t";
                for (my $run=0; $run<=119; $run++)
                {
                    print {$out_pattern[$out_strand]} "$head[$run]";
                    print {$out_pattern[$out_strand]} "\t" if ($run<119);
                }
                print {$out_pattern[$out_strand]} "\n";
            }
        }        
        print "DONE\n";
    }
}
close $out_pattern[0];
close $out_pattern[1];
}

sub option
{
    print "Transcripts_coverage.pl\nChongyuan Luo, Biotech, Rutgers\nSep. 2009\n\nOptions\n";
    print "-f Coverage map corresponding with Watson strand\n";
    print "-r Coverage map corresponding with Crick strand\n";
    print "-d output path\n";
}
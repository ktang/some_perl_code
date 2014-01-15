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

open cdt_file, "$m" or die $!;
my @outputname=split(/cdt/,$m);
my $outputname=$outputname[0]."constance";
open plain_file, ">$outputname" or die $!;

$_=<cdt_file>; chop $_;
my @read;
@read=split(/\t/,$_);
for (my $run=0; $run<=@read-4; $run++)
{
    print plain_file "$run";
    print plain_file "\t" if ($run<@read-4);
}
print plain_file "\n";
$_=<cdt_file>;

while (<cdt_file>)
{
    chop $_;
    @read=split(/\t/,$_);
    print plain_file "$read[1]\t";
    for (my $run=4; $run<=@read-1; $run++)
    {
        print plain_file "$read[$run]";
        print plain_file "\t" if ($run<@read-1);    
    }
    print plain_file "\n";
}
close cdt_file;
close plain_file;

}
sub option
{
    print "cdt_to_plain.pl\nChongyuan Luo, Biotech, Rutgers\nOct. 2009\n\nOptions\n";
    print "-m input pattern file\n";
    print "-d output path\n";
}


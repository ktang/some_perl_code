#!/usr/bin/perl -w

use strict;

#START HERE MAIN PROGRAM

my $usage = "$0 <Ecker_data> <peaks_gff>";
die $usage unless (@ARGV == 2);

$ARGV[0] =~ /(\S+)\.gff/;
my $pre_e = $1;
my $out_e = join ("",$pre_e,"_sorted.gff");

$ARGV[1] =~  /(\S+)\.gff/;
my $pre_p = $1;
my $out_p =  join ("",$pre_p,"_sorted.gff"); 

system ("sortgfffilev03.sh $ARGV[0] $out_e");
#system("killgffcomments.sh sorted.temp $out_e");
system ("sortgfffilev03.sh $ARGV[1] $out_p");
#system("killgffcomments.sh sorted.temp $out_p");
#system("rm -f sorted.temp");
my $ecker_file = $out_e;
my $peakfilenm=$out_p;

system("perl find_num_of_overlap_region_v4.pl $out_e $out_p");

print STDERR  "\a";
exit;
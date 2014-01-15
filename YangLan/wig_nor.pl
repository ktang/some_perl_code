#!usr/bin/perl -w
use strict;
die "Usage:<sam.uniq> <wig> <approximate libsize>\n" unless @ARGV==3;

my $sum=$ARGV[0];
my $lib=$ARGV[2];
=hide
open IN,$ARGV[0];
  while(<IN>){
    next if(/^@/);
    $sum++;
  }
close IN;
=cut

open IN1,$ARGV[1];
  while(<IN1>){
    chomp;
	my @c=split /\t/,$_;
	if(/^track/){
	   print "$_\n";
	}else{
	   my $nor=$c[3]*$lib*1000000/$sum;
	   $c[3]=$nor;
	   my $out=join "\t",@c;
	   print "$out\n";
	}
  }
close IN1;

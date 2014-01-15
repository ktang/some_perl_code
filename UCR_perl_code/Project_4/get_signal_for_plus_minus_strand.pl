#!/usr/bin/perl -w

=head2

this script use 

/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/33_probes_mapped_to_minus_stands.txt
and
/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/data/gff_ten

as input get signal for "-" probe.
=cut
use strict;

my $usage = "$0 <signal_gff>";

die $usage unless (@ARGV == 1);


my $minus_file = "/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/33_probes_mapped_to_minus_stands.txt";

my $signal_file = $ARGV[0];
open (MINUS,  "<$minus_file")
  or die "Cannot open $minus_file: $!";
  
open (IN, "<$signal_file" )
	 or die "Cannot open $signal_file: $!";
  
 my $line;
 
 my %strands;
 
 while ($line = <MINUS>)
 {
 	chomp $line;
 	my @pts = split "\t", $line;
 	$strands{$pts[0]} ->{$pts[1]} =  $pts[3];
 }
  
  close(MINUS);
  
   while ($line = <IN>)
 {
 	chomp $line;
 	my @pts = split "\t", $line;
 	if (exists $strands{$pts[0]}->{$pts[3]})
 	{
 		my $out = $strands{$pts[0]}->{$pts[3]};
 		print "$pts[0]\t$pts[3]\t$pts[4]\t$out\t$pts[5]\n";
 	}
 	
 	
 }
  close(IN);
exit;


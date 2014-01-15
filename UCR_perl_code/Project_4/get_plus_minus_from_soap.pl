#!/usr/bin/perl -w

=head2
finished
this script use 

/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/sorted_position_ID_v2
and

/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/Arab_thal_H1_all_meth_probes_vs_Tair9_chr_all_uniq_exactly_match.soapout

as input get +/- for each probe.
=cut
use strict;

my $soap_file = "/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/Arab_thal_H1_all_meth_probes_vs_Tair9_chr_all_uniq_exactly_match.soapout";

my $probe_file = "/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/sorted_position_ID_v2";
open (SOAP,  "<$soap_file")
  or die "Cannot open $soap_file: $!";
  
open (IN, "<$probe_file" )
	 or die "Cannot open $probe_file: $!";
  
 my $line;
 
 my %strands;
 
 while ($line = <SOAP>)
 {
 	chomp $line;
 	my @pts = split "\t", $line;
 	$strands{$pts[0]} = $pts[6];
 }
  
  close(SOAP);
  
   while ($line = <IN>)
 {
 	chomp $line;
 	my @pts = split "\t", $line;
 	if (exists $strands{$pts[3]})
 	{
 		$pts[4] = 	 $strands{$pts[3]};
 		my $out_line = join "\t",@pts;
 		print "$out_line\n";
 	}
 	
 	else
 	{
 		print STDERR "$line\n";
 		die;
 	}
 }
  close(IN);
exit;


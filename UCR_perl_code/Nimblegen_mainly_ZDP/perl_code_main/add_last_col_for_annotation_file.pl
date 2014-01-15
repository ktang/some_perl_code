#!/usr/bin/perl -w

use strict;
=head2
because the annotated file lack the last column Intergenic.
this scipt use the origanl annotated file to add the Intergenic
col in the new file.
=cut


#START HERE MAIN PROGRAM
#debug this is the debug flag! if it is 1, then debug output

my $usage = "$0 <orignal_annotated file> <lack_last_col_file>";
die $usage unless (@ARGV == 2);

my $orignal_file = $ARGV[0];
my $lack_file = $ARGV[1];

open (ORI, "$orignal_file")
	or die "cannot open $orignal_file:$!";
	
open (LACK, "$lack_file")
	or die "cannot open $lack_file:$!";
	
my %hs;

while (my $line = <ORI>)
{
	chomp $line;
	my @pts = split "\t", $line;
	$hs{$pts[0]}->{$pts[1]} =$pts[11];
		
}
close(ORI);
while (my $line = <LACK>)
{
	chomp $line;
	my @pts = split "\t", $line;
	if (exists $hs{$pts[0]}->{$pts[1]})
	{
		$pts[12] = $hs{$pts[0]}->{$pts[1]};
		my $out = join "\t",@pts;
		print $out,"\n";
	}
	
	else 
	{
		print STDERR "die $line\n\n";
		die;	
	}
}
close (LACK);
exit;
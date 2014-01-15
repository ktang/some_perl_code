#!/usr/bin/perl -w

=head
 not finished!!
=cut

use strict;

sub isoverlappingslop
{
	my $extrabit=$_[0]; #allows for sloppiness - if this is set to zero, then there is no sloppiness
	my $acor=$_[1]-$extrabit;
	my $bcor=$_[2]+$extrabit;
	my $ccor=$_[3]-$extrabit;
	my $dcor=$_[4]+$extrabit;
	
	if ($acor>=$ccor && $acor<=$dcor)
	{return 1;}
	if ($bcor>=$ccor && $bcor<=$dcor)
	{return 1;}
	if ($ccor>=$acor && $ccor<=$bcor)
	{return 1;}
	if ($dcor>=$acor && $dcor<=$bcor)
	{return 1;}
	return 0;
}

#short for bed, long for gff
my $usage = "$0 <bed_file> <anno_gff_file> <out>";
die $usage unless(@ARGV == 3);

my $bed_file=$ARGV[0]; #input filenam
my $anno_gff_file=$ARGV[1]; #second filename
my $output = $ARGV[2]; #extra bit to add to each piece to allow near overlaps

my @gff;
my @pts_bed;
my $line;
my ($gstart,$gend,$bstart,$bend);

print "$long_num,$short_num,$overlap\n";

open (BED,"<$bed_file")
	or die "cannot open file $bed_file:$!";

open (GFF,"<$anno_gff_file")
	or die "cannot opne file $anno_gff_file:$!";

open (OUT,">$output")
	or die "cannot open file $output:$!";

my @gff=<GFF>;
close(GFF);

while ($line = <BED>)
{		
		chomp $line;
		@pts_bed = split "\t", $line;
		
	LOOKLOOP: for ($j=0; $j <= $#gff; $j++)
	{
		my @pts_gff = split "\t", $gff[$j];
		if ( (lc($pts_bed[0]) eq lc($pts_gff[0])) and ($pts_gff[2] ne "chromosome"))
		{
			$gstart=$pts_gff[3];
			$gend=$pts_gff[4];
			$bstart = $pts_bed[1];
			$bend = $pts_bed[2];
			if (isoverlappingslop(0, $gstart,$gend,$bstart,$bend))
			{
				if ( $pts_gff[2] eq "gene")
				{next LOOKLOOP;}
				
				last LOOKLOOP;
			}
		}
	}
}
close(LONG);

my $gaps = 2 * $slopval;

print "long file contains $long_num lines \n short files contains $short_num lines\n overlap within $gaps bp is $overlap\n\a\n";
exit;

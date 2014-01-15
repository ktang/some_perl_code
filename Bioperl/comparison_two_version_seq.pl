#!/usr/bin/perl -w
# use bioperl to get the 
#length of each transcript

use strict;
use Bio::SeqIO;

my $usage = "$0 <our_version> <rebaca_version>";

die $usage unless(@ARGV == 2);

my %genes;
my $name;
my ($same,$diff) = (0,0);
my $old = new Bio::SeqIO(-format => 'fasta',
             		-file => $ARGV[1]);

my $new = new Bio::SeqIO (-format => 'fasta',
			-file => $ARGV[0]);

while( my $oldin = $old->next_seq() )
{
  $name = $oldin->id;
  $genes{$name}= $oldin->seq;
}

$old->close;

while (my $newin = $new->next_seq())
{
	$name = $newin->id;
	if(exists $genes{$name})
	{
		if ($newin->seq eq $genes{$name})
		{
			$same ++;
		}

		else
		{
			$diff++;
			print "$name\n";
		}
	}

}

print "same=$same\ndiff=$diff\n";
print"\a";
exit;

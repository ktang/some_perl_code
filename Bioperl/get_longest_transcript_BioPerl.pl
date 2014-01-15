#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 <longest_trans_names> <cdna> <output>";

die $usage unless(@ARGV == 3);


my %longest;
my $name;

open (IN,"<$ARGV[0]")
	or die "cannot open the file:$!";
	
	
while (my $line = <IN>)
{
	chomp $line;
	$longest{$line} = -1;	
}

 my $seq_in = Bio::SeqIO->new('-file' => "<$ARGV[1]",
                             '-format' => "fasta");
         
  my $seq_out = Bio::SeqIO->new('-file' => ">$ARGV[2]",
                               '-format' => 'fasta');


while( my $seq = $seq_in->next_seq() )
{
	if(exists $longest{$seq->id})
	{
		$seq_out->write_seq($seq)
	}
}

$seq_in->close;
$seq_out->close;
print "\a";
exit;
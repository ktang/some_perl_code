#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 <in1> <out>";
die $usage unless(@ARGV >= 2);
my ($input, $outFile) = @ARGV[0..1];

my $IN1;

my $OUT;


   


open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";


my %gene_len;
my $name;
my $altsp;



my $in = new Bio::SeqIO(-format => 'fasta',
             		-file => $input);
while( my $seq = $in->next_seq() )
{
  
  $name = $seq->id;
  $gene_len{$name}= $seq->length;
 }


foreach my $gene(sort keys %gene_len)
{
  print $OUT $gene,"\t",$gene_len{$gene},"\n";
}

## printf $FNAFILE "Query= $result->query_name  Hit=$result->num_hits \n";
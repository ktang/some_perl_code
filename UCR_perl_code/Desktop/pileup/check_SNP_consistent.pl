#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 <in1> <in2> <out>";
die $usage unless(@ARGV >= 3);
my ($input,$in, $outFile) = @ARGV[0..2];

my $IN;

my $OUT;


open ($IN, '<', $in)
	or die "Cannot open $in:$!";


open ($OUT, '>', $outFile)
  or die "Cannot open $outFile: $!";
  


my $line;
my %genes;
my %dna;
my $bp;
my @parts;
my ($gene_count, $transcript_count, $incon_gene, $incon_transcript)
    = (0,0,0,0);
my %gen;

while ($line = <$IN>)
{
	chomp($line);
	@parts = split /\t/,$line;
	$genes{$parts[0]}->{$parts[1]} = $parts[7];
	$dna{$parts[0]}->{$parts[1]} = $parts[2];
}



my $read = new Bio::SeqIO(-format => 'fasta',
             		-file => $input);
while( my $seq = $read->next_seq() )
{
  if(exists $genes{$seq->id})
  {
  	$transcript_count++;
  	
  	my $gene_name = substr($seq->id,0,9);
     	 
     	 if (exists $gen{$gene_name} )
     	 {
     	    ;
     	 }
     	 else
     	 {
     	 	$gen{$gene_name} = -1;
     	 	$gene_count++;
     	 }
  	
  	
  	
  	foreach my $position(sort keys %{$genes{$seq->id}})
  	{
 	 $bp = substr($seq->seq(),$position-1,1);
 	 if ($bp eq $genes{$seq->id}->{$position})
 	 {
 	    
 	    
 	 
 	 	print $OUT $seq->id,"\t",$position,"\t",
 	 	$dna{$seq->id}->{$position}."\t",
 	 	$genes{$seq->id}->{$position},
 	 	"\t",$bp,"\n";
 	 }
 	}
  }
 }
 
 print "total transcripts in C24_lib:",$transcript_count,"\n";
 print "total genes in C24_lib:",$gene_count,"\n";
 
 
 print "\a";
 exit;
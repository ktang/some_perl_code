#!/usr/bin/perl -w

#this script use Bioperl to read chromosome seq in,
# use tr to count CG number in window size wsize.

use strict;
use Bio::SeqIO;

my $usage = "$0 <chr_fasta> <output_file> <window_size> <chrX>";

die $usage unless(@ARGV == 3);

my $chr_in = Bio::SeqIO -> new ('-file' => "<$ARGV[0]",
			    '-format' => "fasta");
open (OUT,  ">$ARGV[1]")
  or die "Cannot open $ARGV[1]: $!";

my $chr = $ARGV[3]; 
my $chr_seq;
my @parts;
my $line;
my %strand;
my $wsize = $ARGV[2];

$chr_seq_obj = $chr_in->next_seq;
$seq = $chr_seq_obj -> seq;
my $seq_length = length($seq);
my $num_window = int($seq_length/$wsize) + 1;

for (my $i = 1; $i < $num_window; $i++)
{
	my $win_seq = substr($seq, 0, $wsize);
	substr($seq, 0, $wsize) = "";
	my $CG_count = ($win_seq =~ s/CG/CG/g);
	my $start = $wsize*($i-1) +1;
	my $end = $wsize*$i;
	print  OUT "$chr\t.\t.\t$start\t$end\t$CG_count\t.\t.\t.\n";
}
  
#the last window
	my $CG_count = ($seq =~ s/CG/CG/g);
	my $start = $wsize*($num_window - 1) +1;
	my $end = $seq_length;
	print  OUT "$chr\t.\t.\t$start\t$end\t$CG_count\t.\t.\t.\n";

  close (OUT);
  print STDERR "\a";
  exit;

#!/usr/bin/perl -w

#this script use Bioperl to read chromosome seq in,
# use tr to count CG number in window size wsize, 
# but the window moves $move_size forward.

use strict;
use Bio::SeqIO;

my $usage = "$0 <chr_fasta> <output_file> <window_size> <move_seze>";

die $usage unless(@ARGV == 4);

my $chr_in = Bio::SeqIO -> new ('-file' => "<$ARGV[0]",
			    '-format' => "fasta");
open (OUT,  ">$ARGV[1]")
  or die "Cannot open $ARGV[1]: $!";

my $wsize = $ARGV[2];
my $move_size = $ARGV[3];

while (my $chr_seq_obj = $chr_in->next_seq)
{
	my $chr_name = lc($chr_seq_obj -> id);

	if ($chr_name eq "chrc")
	{	
		$chr_name = "chrC";
	}
	elsif ($chr_name eq "chrm")
	{
		$chr_name = "chrM";
	}
	
	my $chr_seq = $chr_seq_obj -> seq;
	my $seq_length = length($chr_seq);
	my $num_window = int($seq_length/$wsize) + 1;
	my ($start,$end, $i);
	for ($i = 0; ($wsize + $i * $move_size) <= $seq_length; $i++)
	{
		$start = $i * $move_size +1;
		$end = $start + $wsize - 1;
		my $win_seq = substr ($chr_seq, $start - 1, $wsize);
		my $CG_count = ($win_seq =~ s/CG/CG/g);
		if ($CG_count = 0)
		{
			$CG_count = 0;
		}
		my $last = substr($winz_seq,499,1);
		my $first = substr($chr_seq , $end,1);
		if ( ($last eq "C") and ($first eq "G"))
		{
			$CG_count++;
		}
		print  OUT "$chr\t.\t.\t$start\t$end\t$CG_count\t.\t.\t.\n";
	}

	#maybe the last window
	$start = $i * $move_size + 1;
	$end = $seq_length;
	$len = $end - $start + 1;
	if ( $len >= $wsize)
	{
		print "ERROR:start =$start, end =$end, i = $i\n";
	} 
	my $win_seq = ($chr_seq, $start - 1, $len)
	my $CG_count = ($win_seq =~ s/CG/CG/g);
	if ($CG_count == 0)
	{
		$CG_count = 0;
	}	
	print  OUT "$chr\t.\t.\t$start\t$end\t$CG_count\t.\t.\t.\n";
}

close (OUT);
print STDERR "\a";
exit;

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
my ($CG_num,$CHG_num,$CHH_num) = (0,0,0);

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
	 	($CG_num,$CHG_num,$CHH_num) = (0,0,0);
		$start = $i * $move_size +1;
		$end = $start + $wsize - 1;
		my $win_seq = substr ($chr_seq, $start - 1, $wsize+2);
		my @win_seq_array = split "",$win_seq;
		
		for (my $j = 0; $j <= 499; $j++)
		{
			if ($win_seq_array[$j] eq "C")
			{
				if ($win_seq_array[$j+1] eq "G" )
				{$CG_num++}
				else
				{
					if($win_seq_array[$j+2] eq "G")
					{
						$CHG_num++;
					}
					else{$CHH_num++}
				}
			}
		} 
		print  OUT "$chr_name\t.\t.\t$start\t$end\t$CG_num\t$CHG_num\t$CHH_num\t.\n";
	}

	#maybe the last window
	
	($CG_num,$CHG_num,$CHH_num) = (0,0,0);
	$start = $i * $move_size + 1;
	$end = $seq_length;
	my $len = $end - $start + 1;
	if ( $len >= $wsize)
	{
		print "ERROR:start =$start, end =$end, i = $i\n";
	} 
	my $win_seq = substr($chr_seq, $start - 1, $len);
	
	my @win_seq_array = split "",$win_seq;

	for (my $j = 0; $j < $len-2 ; $j++)
	{
		if ($win_seq_array[$j] eq "C")
		{
			if ($win_seq_array[$j+1] eq "G" )
			{$CG_num++}
			else
			{
				if($win_seq_array[$j+2] eq "G")
				{
					$CHG_num++;
				}
				else{$CHH_num++}
			}
		}
	}
	if ( ($win_seq_array[$len-2] eq "C") and ($win_seq_array[$len-1] eq "G"))
	{$CG_num++} 
	print  OUT "$chr_name\t.\t.\t$start\t$end\t$CG_num\t$CHG_num\t$CHH_num\t.\n";
}

	
close (OUT);
print STDERR "\a";
exit;

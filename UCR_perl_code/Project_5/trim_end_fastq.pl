#! /usr/bin/perl -w
#not finished, FASTX can do it.
=head2
Copyright (C) 2010 Kai Tang

Nov 22, 2010
Noet: in variant = "illumina" B = 2;
in variatn = "sanger" which is default, # is 2, B = 33,T = 20
my $debug = 1;
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Seq::Quality;
my $usage = "$0 <input_fastq> <output_fastq> <stat_file> <variant>";

my $min_len = 20;
die $usage unless (@ARGV == 4);
my ( $in_fq, $out_fq, $stat_file, $variant) = @ARGV[0..3]; 
my $seq_in = Bio::SeqIO->new(-format=>'fastq',
			      -variant=> $variant,
				-file=>$in_fq);

my $seq_out = Bio::SeqIO->new(-format=>'fastq',
			       -file=>">$out_fq" );

open ( OUT , ">$stat_file")
	or die "cannot open stat_file $stat_file:$!";

my $no_total = 0;
my %trims;
my $no_not_trim = 0;
while (my $inseq = $seq_in->next_seq)
{
	$no_total ++;
#	if ($debug) { print STDERR "qual:",$inseq->qual_text,"\n";}
	my $score = $inseq->qual_text;
	my @pts = split "  ",$score;
	my $score2 = join"",@pts;
#	if ($debug)
#	{print STDERR "score:$score\n";}
	my $len = length($score);
	if ($debug) {my $l = length($score2);print STDERR "\n\nlength:$l \n$score2\n";}
	my $left_len= $len;
	for (my $i = $len - 1; $i >= $min_len - 1;$i--)
	{
		my $thisScore = substr($score,$i,1);
		if ($thisScore < 20){
			$left_len--;
		}		
		else{
			last;
		}
	}
	if ($left_len< $min_len)
	{}
	else{
		my $new_seq = substr($inseq->seq,0,$left_len);
		my $new_score = substr($score,0,$left_len);
		my $seq_fq = Bio::Seq::Quality->new(-id=>inseq->id,-seq=>$new_seq,-qual=>\$new_score);
		$seq_out->write_fastq($seq_fq);
	}	
}
=cut

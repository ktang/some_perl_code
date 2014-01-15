#!/usr/bin/perl -w

# use method in PNAS 97:11911
# 80 consecutive residues containing at least 30 Q and/or N

#input /Users/tang58/DataBase/TAIR10/Proteins/TAIR10_pep_20101214

use strict;
use File::Spec;
use Bio::SeqIO;

my $win_size = 80;
my $cutoff = 30;

my $cutoff_per = $cutoff / $win_size;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

my $infileformat = "fasta";
my $outfileformat = "fasta";

my $seq_in = Bio::SeqIO->new('-file' => "<$input",
                             '-format' => $infileformat);
my $seq_out = Bio::SeqIO->new('-file' => ">$output",
                              '-format' => $outfileformat);



# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
	my $seq = $inseq->seq;
	
	my $length =  length ($seq) ;
	
	if ( $length < $win_size ){
		
		my $sub_seq = $seq;
		my $Q_num = ( $sub_seq =~ tr/Q/Q/);
		my $N_num = ( $sub_seq =~ tr/N/N/);
		if( ($Q_num  + $N_num ) / $length >= $cutoff_per){
			$seq_out->write_seq($inseq);
			print STDERR $inseq->id, "\n";
		}
	}else{
		for my $i(0..($length - $win_size)){
			my $sub_seq = substr ($seq, $i, $win_size);
			my $Q_num = ( $sub_seq =~ tr/Q/Q/);
			my $N_num = ( $sub_seq =~ tr/N/N/);
			if( $Q_num  + $N_num >= $cutoff){
				$seq_out->write_seq($inseq);
				last;
			}
		}
	}
}
exit;
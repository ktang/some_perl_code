#!/usr/bin/perl -w

my $debug = 0;
if($debug){
	print STDERR "debug = 1 \n\n";
}

use File::Spec;
use strict;
use Bio::SeqIO;
my $usage = "$0 <indir> <input> <outdir> <out_prefix>";
die $usage unless (@ARGV == 4);

my $indir = shift or die "indir";
my $input = shift or die "input";
my $outdir = shift or die "outidr";
my $outpre = shift or die "outpre";

my $infile = File::Spec->catfile($indir, $input);

die unless(-e $infile );



#	my $seq_out = Bio::SeqIO->new(-file => ">$outfile",
#								-format => 'fasta',
#								-width => 80);

my $seq_in = Bio::SeqIO->new(-file   => $infile,
							 -format => "fasta");
							 
while(my $seq = $seq_in->next_seq){
	my $chr = lc ($seq->id);
	my $output = File::Spec->catfile($outdir, $outpre . "_". $chr . ".fa");
	die "$output \n" if (-e $output);
	
	if($debug) {
		print STDERR $output, "\n";
	}
	else{
		my $seq_out = Bio::SeqIO->new(-file   => ">$output",
									  -format => 'fasta',
									  -width  => 60);
#		$seq->id = $chr;
		my $outseq = Bio::Seq->new(-seq => $seq->seq(),
									-id  => $chr);
		$seq_out->write_seq($outseq);
	}
}


exit;
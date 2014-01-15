#!/usr/bin/perl -w
# use bioperl to do 
# bl2seq

use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SimpleAlign;

my $usage = "$0 <Col-0_cdna> <C24_cdna>";
die $usage unless(@ARGV == 2);

my %genes_Col0;
my %genes_C24;

my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastn',
                                               'outfile' => "bl2seq.out");

my $str = Bio::AlignIO->new(-file=> ">bl2seq.out",'-format' => 'bl2seq');


my $inseq1 = Bio::SeqIO->new(-file   => "<$ARGV[0]",
                            -format => 'fasta' );
    
while (my $seq1 = $inseq1->next_seq) 
{
    $genes_Col0{$seq1->id} = $seq1;     
}

my $inseq2 = Bio::SeqIO->new(-file   => "<$ARGV[1]",
                            -format => 'fasta' );

while (my $seq2 = $inseq2->next_seq)
{
	if (exists $genes_Col0{$seq2->id} )
	{
		my $bl2seq_report = $factory->bl2seq($genes_Col0{$seq2->id}, $seq2);
		$str->write_aln($bl2seq_report);
	}
}
exit;

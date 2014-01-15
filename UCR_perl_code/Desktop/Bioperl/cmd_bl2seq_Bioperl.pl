#!/usr/bin/perl -w
# use bioperl to do 
# bl2seq

use strict;
use Bio::SeqIO;

#use Bio::SearchIO;
#use Bio::AlignIO;
#use Bio::Tools::Run::StandAloneBlast;
#use Bio::SimpleAlign;

my $usage = "$0 <Col-0_cdna> <C24_cdna>";
die $usage unless(@ARGV == 2);

my %genes_Col0;
my %genes_length;

#my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastn',
#                                               'outfile' => "bl2seq.out");

#my $str = Bio::AlignIO->new(-file=> ">bl2seq.out",'-format' => 'bl2seq');


my $inseq1 = Bio::SeqIO->new(-file   => "<$ARGV[0]",
                            -format => 'fasta' );
    
while (my $seq1 = $inseq1->next_seq) 
{
    $genes_Col0{$seq1->id} = $seq1;
    $genes_length{$seq1->id} = $seq1->length();  
}

$inseq1->close;

my $inseq2 = Bio::SeqIO->new(-file   => "<$ARGV[1]",
                            -format => 'fasta' );

while (my $seq2 = $inseq2->next_seq)
{
	if (exists $genes_Col0{$seq2->id} )
	{
			#open(SUBJECT,">subject_file") or die "cannot open query file";
			#	print SUBJECT $genes_Col0{$seq2->id},"\n";
			#close(SUBJECT);
			
			my $sub_seq = Bio::SeqIO->new(-file   => ">subject_file",
                            -format => 'fasta' );
            $sub_seq->write_seq($genes_Col0{$seq2->id} );
            $sub_seq->close;
			
			my $query_seq = Bio::SeqIO->new(-file   => ">query_file",
                            -format => 'fasta' );
            $query_seq->write_seq($seq2);
            $query_seq->close;
            
            my $len =  $genes_length{$seq2->id};
        my $cmd = "blastn -task blastn -query ./query_file -subject ./subject_file -subject_loc 0-$len -out temp_bl2seq_output -outfmt 6 -evalue 1e-10 -num_alignments 1 -max_target_seqs 1";
        `$cmd`;
        `cat temp_bl2seq_output >> bl2seq_C24_snponly_vs_Col0_output`;
	}
}

$inseq2->close;
exit;
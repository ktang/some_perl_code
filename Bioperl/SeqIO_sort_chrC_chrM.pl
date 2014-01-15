#!/usr/bin/perl -w

use Bio::SeqIO;
use Bio::Seq;
use strict;
         
         # get command-line arguments, or die with a usage statement
         my $infile = $ARGV[0];
         my $outfile = $ARGV[1];
         # create one SeqIO object to read in,and another to write out
         my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                                      '-format' => 'fasta');
            
          my $seq_out ;                            
         if (-e $outfile)
         {	
         		die "file $outfile exists!!!\n";
         		
         }
         
         else 
         {
            $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                                       '-format' => 'fasta');
         }
		 my %hash;
         # write each entry in the input file to the output file
         while (my $inseq = $seq_in->next_seq) {
         	my $id = $inseq->id;
         	my $post = substr($id,3,1);
         	my $key = "Chr".$post;
         	$hash{$key} = $inseq->seq;
            
         }
         
         foreach my $key (sort keys %hash)
         {
         		my $seq = Bio::Seq->new(-id => $key,-seq => $hash{$key} );
         		$seq_out->write_seq($seq);
         }
         
         exit;
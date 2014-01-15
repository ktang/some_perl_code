#!/usr/bin/perl -w

# get unique seqs
use strict;
use File::Spec;
use Bio::SeqIO;

my $usage = "$0 <input_file>";
die $usage unless(@ARGV >= 1);
my $inFile = $ARGV[0];

my $inDir = "/mnt/disk2/kai/analysis_xingyu/xingyu_RNA_seq_FASTA_removed";
my $outDir =  "/mnt/disk2/kai/analysis_xingyu/unique";

	my @parts = split /\./, $inFile;
	my $ext = pop @parts;
	my $outFile = File::Spec->catfile($outDir, join(".", @parts) . "_unique" . "." . $ext);
	my $input = File::Spec->catfile($inDir, $inFile);
    my $seqin = Bio::SeqIO->new(-file=>$input, -format=>'fasta');
    my %seqs;
	my %len;
    my ($min_len, $max_len) = (100000000000, 0);
    my $n = 0;
    while(my $seq = $seqin->next_seq){
	if(defined $seqs{$seq->seq}->{$seq->id}){
		die "sequence ", $seq->id, ": ", $seq->seq, " is redundant";
	}else{
		$n++;
	    $seqs{$seq->seq}->{$seq->id} = $seq;
		if($seq->length < $min_len){
			$min_len = $seq->length;
		}
		if($seq->length > $max_len){
			$max_len = $seq->length;
		}
	}
	}
	my $seqout = Bio::SeqIO->new(-file=>">$outFile", -format=>'fasta');
	foreach my $s(keys %seqs){
		my @seqIds = sort keys %{$seqs{$s}};
		my $repSeq = $seqs{$s}->{$seqIds[0]};
		$repSeq->desc("RDN=" . scalar(@seqIds));
		$seqout->write_seq($repSeq);
		$len{$repSeq->length}++;
	}
	$seqin->close;
	$seqout->close;


    my $out_sum = "$outDir/$parts[0]_summary";
    my $SUM;
open ($SUM, '>', $out_sum)
  or die  "Cannot open:$out_sum:$!";
     print $SUM "In query file $inFile, Number seqs: $n, max length: $max_len, min_len: $min_len, number of unique seqs: ", scalar(keys %seqs), "\n";

    print $SUM "Unique sequence length distribution:\n";
    foreach my $l(sort {$a<=>$b} keys %len){
	    print $SUM $l, "\t", $len{$l}, "\n";
    }
  close($SUM);
    print STDERR "In query file $inFile, Number seqs: $n, max length: $max_len, min_len: $min_len, number of unique seqs: ", scalar(keys %seqs), "\n";

    print STDERR "Unique sequence length distribution:\n";
    foreach my $l(sort {$a<=>$b} keys %len){
	    print STDERR $l, "\t", $len{$l}, "\n";
    }



print "\a";
exit;

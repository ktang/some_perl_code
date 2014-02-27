#!/usr/bin/perl -w
## given a list of mature miRNA seqs, count number of reads in each library
## that had perfect match to each
use strict;
use Bio::SeqIO;

print STDERR "STDOUT \n\n";

# /Users/tang58/DataBase/miRBase/Release_19/ath_mature_r19.fa
my $ten_million = 10000000;
my $usage = "$0\n <mature miRNAs> <format> <clean reads in lib1> [<lib2>.....] \n\n ";
die $usage unless(@ARGV >= 3);
my $miRNA_file = shift or die;
my $format = shift or die;

die "fastq or fasta" unless ( $format eq "fasta" or $format eq "fastq") ; 

my $seqin = Bio::SeqIO->new(-file=>$miRNA_file, -format=>'fasta');
my %miRNAs;
while(my $seq = $seqin->next_seq){
	my $seqstr = $seq->seq;
	$seqstr =~ tr/U/T/;
	if(!defined $miRNAs{$seqstr}){
		$miRNAs{$seqstr} = $seq->id;
	}else{
		$miRNAs{$seqstr} = $miRNAs{$seqstr} . ";" . $seq->id;
	}
}
$seqin->close;

my %miRNA_counts;
my %miRNA_normal; # normalized to TPTM
foreach my $file(@ARGV){
	#my $in = Bio::SeqIO->new(-file=>$file, -format=> $format );
	my $in = Bio::SeqIO->new(-file=>"zcat $file|", -format=> $format );
	my $total = 0;
	while(my $seq = $in->next_seq){
		if(defined $miRNAs{$seq->seq}){
			$miRNA_counts{$file}->{$seq->seq}++;
		}
		$total++;
	}
	foreach my $k(keys %miRNAs){
		if(defined $miRNA_counts{$file}->{$k}){
			$miRNA_normal{$file}->{$k} = $miRNA_counts{$file}->{$k} * $ten_million / $total;
		}else{
			$miRNA_counts{$file}->{$k} = 0;
			$miRNA_normal{$file}->{$k} = 0;
		}
	}
	$in->close;
}
print join("\t", ("miRNAId", "miRNASeq"));
foreach my $file(@ARGV){
	print "\t", "Raw_".$file;
}
foreach my $file(@ARGV){
	print "\t", "TPTM_".$file;
}
print "\n";
foreach my $k(keys %miRNAs){
	print $miRNAs{$k}, "\t", $k;
	foreach my $file(@ARGV){
		print "\t", $miRNA_counts{$file}->{$k};
	}
	foreach my $file(@ARGV){
		print "\t", $miRNA_normal{$file}->{$k};
	}
	print "\n";
}

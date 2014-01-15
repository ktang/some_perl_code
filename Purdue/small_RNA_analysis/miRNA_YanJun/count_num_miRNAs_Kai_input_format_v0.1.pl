#!/usr/bin/perl -w
## given a list of mature miRNA seqs, count number of reads in each library
## that had perfect match to each

#v0.1
# used for shift miRNA
use strict;
use Bio::SeqIO;

# /Users/tang58/DataBase/miRBase/Release_19/ath_mature_r19.fa
my $ten_million = 10000000;
#my $usage = "$0 <mature miRNAs> <format> <clean reads in lib1> [<lib2>.....] ";
#die $usage unless(@ARGV >= 3);
#my $miRNA_file = shift or die;
#my $format = shift or die;
#die "fastq or fasta" unless ( $format eq "fasta" or $format eq "fastq") ; 

#  /Users/tang58/DataBase/miRBase/Release_19/Kai/YanJun_all_possible_shift_seq_db.txt

my $usage = "$0 \n<shfit_miRNAs_list> <format> <clean reads in lib1> [<lib2>.....] \n\n ";
die $usage unless(@ARGV >= 3);

my $shfit_miRNAs_list = shift or die;
my $format = shift or die;
die "fastq or fasta" unless ( $format eq "fasta" or $format eq "fastq") ; 

#my $seqin = Bio::SeqIO->new(-file=>$miRNA_file, -format=>'fasta');
#my %miRNAs;
#while(my $seq = $seqin->next_seq){
#	my $seqstr = $seq->seq;
#	$seqstr =~ tr/U/T/;
#	if(!defined $miRNAs{$seqstr}){
#		$miRNAs{$seqstr} = $seq->id;
#	}else{
#		$miRNAs{$seqstr} = $miRNAs{$seqstr} . ";" . $seq->id;
#	}
#}
#$seqin->close;

my %miRNAs;
open (IN, $shfit_miRNAs_list) or die;
while (<IN>){
	chomp;
	my @a = split "\t";
	#$miRNAs{$seqstr} = $seq->id;
	$miRNAs{$a[0]} =  $a[1];
}
close IN;

my %miRNA_counts;
my %miRNA_normal; # normalized to TPTM
foreach my $file(@ARGV){
	my $in = Bio::SeqIO->new(-file=>$file, -format=> $format );
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
print join("\t", ("shift_miRNAId", "shift_miRNASeq"));
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

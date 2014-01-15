#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 only detect SNP in CDS, output the gene name; 

Nov 23,2010
=cut

use strict;
my $usage = "$0 <input_SNP_list>";
die $usage unless(@ARGV == 1);

my $snp_file = $ARGV[0];
my $gff = "/Users/kaitang/Desktop/TAIR/TAIR9_GFF3_genes.gff";
my %snps_h;
my %gene_h;
my %ref_h;
open (SNP,$snp_file)
	or die "cannot open $snp_file:$!";
open (GFF,$gff)
	or die "cannot open $gff:$!";

while (<SNP>){
	chomp;
	my @pts = split /\t/;
	my $chr = lc $pts[0];
	$ref_h{$chr}->{$pts[1]} = $pts[2];
	$snps_h{$chr}->{$pts[1]}->{1}= $pts[3];
	$snps_h{$chr}->{$pts[1]}->{2}= $pts[5];
}

while (<GFF>){
	chomp;
	my @temp = split /\t/;
	my $chr = lc $temp[0];
	if ($temp[2] eq 'CDS'){
		foreach my $pos (keys %{$ref_h{$chr}}){
			if ($pos >= $temp[3] and $pos <= $temp[4]){
				my @parents = split /,/,$temp[8];
				my @genes = split /=/, $parents[0];
				$gene_h{$chr}->{$pos} = $genes[1];
			}
		}
	}
}

foreach my $chr (sort keys %gene_h){
	foreach my $pos (sort {$a <=> $b} keys %{$gene_h{$chr}}){
		my $ref_base = $ref_h{$chr}->{$pos};
		my $best = $snps_h{$chr}->{$pos}->{1};
		my $second =  $snps_h{$chr}->{$pos}->{2};
		my $gene = $gene_h{$chr}->{$pos};
		my $line = join "\t",($chr,$pos,$ref_base,$best,$second,$gene);
		print $line,"\n";
	}
}
exit;

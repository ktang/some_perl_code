#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 2 this time add a function to 
determine which position it is in codon(1,2 or 3) 

Nov 23,2010
=cut

use strict;
my $usage = "$0 <input_SNP_list>";
die $usage unless(@ARGV == 1);

my $snp_file = $ARGV[0];
my $gff = "/Users/kaitang/Desktop/TAIR/TAIR9_GFF3_genes.gff";
my %snps_h;
my %gene_h;
my %ref_h;#store the best and second best snp base;
my %codon_h;#which position in a codon,1,2 or 3
my %counts_h;#store the total count{0},best_count{1};second_best{2}
open (SNP,$snp_file)
	or die "cannot open $snp_file:$!";
open (GFF,$gff)
	or die "cannot open $gff:$!";

while (<SNP>){
	chomp;
	my @pts = split /\t/;
	my $chr = lc $pts[0];
#	if ($pts[2] ne $pts[5] or $pts[12] >= 5){
		$ref_h{$chr}->{$pts[1]} = $pts[2];
		$snps_h{$chr}->{$pts[1]}->{1}= $pts[5];
		$snps_h{$chr}->{$pts[1]}->{2}= $pts[9];
		$counts_h{$chr}->{$pts[1]}-> {0} = $pts[13];
		$counts_h{$chr}->{$pts[1]}-> {1} = $pts[8];
		$counts_h{$chr}->{$pts[1]}-> {2} = $pts[12];
#	}
}

while (<GFF>){
	chomp;
	my @temp = split /\t/;
	my $chr = lc $temp[0];
	if ($temp[2] eq 'CDS'){
		FOR:	foreach my $pos (sort {$a <=> $b} keys %{$ref_h{$chr}}){
			if ($pos >= $temp[3] and $pos <= $temp[4]){
				my @parents = split /,/,$temp[8];
				my @genes = split /=/, $parents[0];
				$gene_h{$chr}->{$pos} = $genes[1];
				if ($temp[6] eq "+"){
					$codon_h{$chr}->{$pos} = ($pos - $temp[3] + 3 - $temp[7]) % 3 + 1;
				}
				elsif ($temp[6] eq "-"){
					$codon_h{$chr}->{$pos} = ($temp[4] - $pos + 3 - $temp[7]) % 3 + 1;
				}
				else{die "neigth + nor -"}
			}
		#	elsif ($pos > $temp[4]){
		#		last FOR;
		#	}
		}
	}
}

foreach my $chr (sort keys %gene_h){
	foreach my $pos (sort {$a <=> $b} keys %{$gene_h{$chr}}){
		my $ref_base = $ref_h{$chr}->{$pos};
		my $best = $snps_h{$chr}->{$pos}->{1};
		my $second =  $snps_h{$chr}->{$pos}->{2};
		my $gene = $gene_h{$chr}->{$pos};
		my $codon_pos = $codon_h{$chr}->{$pos};
		my $total = $counts_h{$chr}->{$pos}-> {0} ;
		my $best_count = $counts_h{$chr}->{$pos}-> {1};
		my $second_count = $counts_h{$chr}->{$pos}->{2};
		my $line = join "\t",($chr,$pos,$ref_base,$total,$best,$best_count,$second,$second_count,$gene,$codon_pos);
		print $line,"\n";
	}
}
exit;

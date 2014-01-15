#!/usr/bin/perl -w
# generate position infor for the position table;
use strict;
use Bio::SeqIO;

my $ref = "/mnt/disk4/renyi/bis_seq/col_ros2_11_09_10/references_names.txt";
my $genome = "Ath_Col0";
my $pos_id = 1;
if(@ARGV >= 3){
	($ref, $genome, $pos_id) = @ARGV[0..2];
}
print join("\t", ("pos_id", "genome", "chr", "pos", "strand", "type")), "\n";
open(RF, $ref) or die "can't open $ref: $!";
while(<RF>){
	chomp;
	if(/\w+/){
		my $file = $_;
		my $seqin = Bio::SeqIO->new(-file=>$file, -format=>'fasta');
		print STDERR "processing $file\n";
		while(my $seq = $seqin->next_seq){
			my $chr = $seq->id;
			my $seq_str = $seq->seq;
			my $chr_len = length($seq_str);
			foreach my $i(0..($chr_len-1)){
				my $pos = $i+1;
				my $base = substr($seq_str, $i, 1);
				my ($strand, $type) = ("+", "CG");
				if($base eq 'C'){
					if(($i+1) <= ($chr_len -1)){
						my $base2 = substr($seq_str, ($i+1), 1);
						if($base2 ne 'G'){
							$type = 'CHG';
							if(($i+2) <= ($chr_len - 1)){
								my $base3 = substr($seq_str, ($i+2),1);
								if($base3 ne 'G'){
									$type = 'CHH';
								}
							}
						}
					}
						print join("\t", ($pos_id, $genome, $chr, $pos, 
						              $strand, $type)), "\n";
						$pos_id++;
				}elsif($base eq 'G'){
					$strand = '-';
					if(($i-1) >= 0){
						my $base2 = substr($seq_str, ($i-1), 1);
						if($base2 ne 'C'){
							$type = 'CHG';
							if(($i-2) >= 0){
								my $base3 = substr($seq_str, ($i-2), 1);
								if($base3 ne 'C'){
									$type = 'CHH';
								}
							}
						}
					}
					print join("\t", ($pos_id, $genome, $chr, $pos, 
						              $strand, $type)), "\n";
						$pos_id++;
				}
			}
		}
	}
}
					


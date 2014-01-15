#!/usr/bin/perl -w
# parse BSMAP output and find methylated positions
use strict;
use Bio::SeqIO;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478);

my %meth;
my %cover;
while(<>){
        chomp;
    my ($id, $seq, $map_flag, $chr, $start, $strand, $num_mm, $mm_info, $mC_loc)        = split /\t/;
        my $str = "+";
        if($strand eq "-+" || $strand eq "--"){
                $str = "-";
        }
        my @temp = split /\:/, $mC_loc;
        if(@temp > 1){
                foreach my $i(1..$#temp){
                        $meth{$chr}->{$strand}->{$temp[$i]}++;
                }
        }
        foreach my $i($start..($start+length($seq)-1)){
                if($strand eq "++" || $strand eq "+-"){
                    $cover{$chr . "+"}->[$i]++;
                }else{
                        $cover{$chr . "-"}->[$i]++;
                }
        }

}

my $chr_file = "/home/renyi/data/ath/tair9/TAIR9_chr_all.fas";
my $seqin = Bio::SeqIO->new(-file=>$chr_file, -format=>'fasta');
my %chrs;
while(my $seq = $seqin->next_seq){
        my @chars = split ("", $seq->seq);
        $chrs{$seq->id} = \@chars;
}

my %num_sites;
my %num_meth_sites;
# do positive strand first
foreach my $chr(sort keys %chr_len){
        foreach my $i(1..$chr_len{$chr}){
                if(defined $cover{$chr . "+"}->{$i} && $cover{$chr . "+"}->{$i} >= 2){
                        if($chrs{$chr}->[$i - 1] ne 'C'){
                                next;
                        }
                        my $type = "CG";
                        if($i < $chr_len{$chr} && $chrs{$chr}->[$i] ne 'G'){
                                $type = "CHG";
                        }
                        if(($i < $chr_len{$chr} - 1) && $chrs{$chr}->[$i+1] ne 'G'){
                                $type = "CHH";
                        }
                        $num_sites{$chr . "+"}->{$type}++;
                        if(defined $meth{$chr}->{"+"}->{$i} && 
                                   $meth{$chr}->{"+"}->{$i} >= 1){
                                print join("\t", ($chr, $i, "+", $type, $meth{$chr}->{"+"}->{$i}, $cover{$chr . "+"}->{$i})), "\n";
                                if($meth{$chr}->{"+"}->{$i} >= 0.5 * $cover{$chr . "+"}->{$i}){
                                        $num_meth_sites{$chr . "+"}->{$type}++;
                                }
                        }
                }
        }
}

# negative strand
foreach my $chr(sort keys %chr_len){
        foreach my $i(1..$chr_len{$chr}){
                if(defined $cover{$chr . "-"}->{$i} && $cover{$chr . "-"}->{$i} >= 2){
                        if($chrs{$chr}->[$i - 1] ne 'G'){
                                next;
                        }
                        my $type = "CG";
                        if($i > 1 && $chrs{$chr}->[$i-2] ne 'C'){
}
                        if(($i > 2) && $chrs{$chr}->[$i-3] ne 'C'){
                                $type = "CHH";
                        }
                        $num_sites{$chr . "-"}->{$type}++;
                        if(defined $meth{$chr}->{"-"}->{$i} && 
                                   $meth{$chr}->{"-"}->{$i} >= 1){
                                print join("\t", ($chr, $i, "-", $type, $meth{$chr}->{"-"}->{$i}, $cover{$chr . "-"}->{$i})), "\n";
                                if($meth{$chr}->{"-"}->{$i} >= 0.5 * $cover{$chr . "-"}->{$i}){
                                        $num_meth_sites{$chr . "-"}->{$type}++;
                                }
                        }
                }
        }
}

## print statistics
print STDERR "Considering only sites with >=2 coverage.\n";
print STDERR "Strand\tChromosome\tType\tNum_C_sites\tNum_mC_sites\tPercentage\n";
foreach my $str("+", "-"){
foreach my $type("CG", "CHG", "CHH"){
        my ($total_C, $total_mC) = (0,0);
    foreach my $chr(sort keys %chr_len){
                $total_C += $num_sites{$chr . $str}->{$type};
                $total_mC += $num_meth_sites{$chr . $str}->{$type};
                print STDERR join("\t", ($str, $chr, $type, $num_sites{$chr . $str}->{$type}, $num_meth_sites{$chr . $str}->{$type}, $num_meth_sites{$chr . $str}->{$type}/$num_sites{$chr . $str}->{$type})), "\n";
        }
}
}

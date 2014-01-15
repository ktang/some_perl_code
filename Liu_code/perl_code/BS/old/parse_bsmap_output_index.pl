#!/usr/bin/perl -w
# parse BSMAP output and find methylated positions
use strict;
use Bio::SeqIO;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478);
my $debug = 0;

my %meth;
my %cover;
while(<>){
        chomp;
    my ($id, $seq, $map_flag, $chr, $start, $strand, $num_mm, $mm_info, $mC_loc)        = split /\t/;
        #next if($map_flag ne 'UM');
        my $str = "+";
        if($strand eq "-+" || $strand eq "--"){
                $str = "-";
        }
        my @temp = split /\:/, $mC_loc;
        if(@temp > 1){
                foreach my $i(1..$#temp){
                        if($debug){
                                #print STDERR "Methylated at $chr, $str, $temp[$i]\n";
                        }
                            $meth{$chr}->{$str}->{$temp[$i]}++;
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
        #$chrs{$seq->id} = \@chars;
        $chrs{$seq->id} = $seq->seq;
}
my $meth_percent_cutoff = 0.9;
my %num_sites;
my %num_meth_sites;
# do positive strand first
foreach my $chr(sort keys %chr_len){
        my $pos = -1;
        while(($pos = index($chrs{$chr}, 'C', $pos)) > -1){
                my $i = $pos + 1;
                if(defined $cover{$chr . "+"}->[$i] && $cover{$chr . "+"}->[$i] >= 2){
                        my $type = "CG";
                        if($i < $chr_len{$chr} - 1){
                                my $triplet = substr($chrs{$chr}, $pos, 3);
                                if($triplet !~ /^CG/){
                                        if($triplet =~ /C\wG/){
                                                $type = "CHG";
                                        }else{
                                                $type = "CHH";
                                        }
                                }
                        }elsif($i < $chr_len{$chr}){
                                my $last_two = substr($chrs{$chr}, -1, 2);
                                if($last_two ne "CG"){
                                        $type = "CHG";
                                }
                        }
                                
                        $num_sites{$chr . "+"}->{$type}++;
                        if(defined $meth{$chr}->{"+"}->{$i} && 
                                   $meth{$chr}->{"+"}->{$i} >= 1){
                                print join("\t", ($chr, $i, "+", $type, $meth{$chr}->{"+"}->{$i}, $cover{$chr . "+"}->[$i])), "\n";
                                if($meth{$chr}->{"+"}->{$i} >= $meth_percent_cutoff * $cover{$chr . "+"}->[$i]){
                                        $num_meth_sites{$chr . "+"}->{$type}++;
}
                        }
                }
                $pos++;
        }
}

# negative strand
foreach my $chr(sort keys %chr_len){
        my $pos = -1;
        while(($pos = index($chrs{$chr}, 'G', $pos)) > -1){
                my $i = $pos + 1;
                if(defined $cover{$chr . "-"}->[$i] && $cover{$chr . "-"}->[$i] >= 2){
                        my $type = "CG";
                        if($i >= 3){
                                my $triplet = substr($chrs{$chr}, $pos-2, 3);
                                if($triplet !~ /CG$/){
                                        if($triplet =~ /C\wG/){
                                                $type = "CHG";
                                        }else{
                                                $type = "CHH";
                                        }
                                }
                if($debug){
                                        print STDERR "Triplet at $i is $triplet, type is $type\n";
                                }

                        }elsif($i == 2){
                                my $first_two = substr($chrs{$chr}, 0, 2);
                                if($first_two ne "CG"){
                                        $type = "CHG";
                                }
                        }
                        $num_sites{$chr . "-"}->{$type}++;
                        if(defined $meth{$chr}->{"-"}->{$i} && 
                                   $meth{$chr}->{"-"}->{$i} >= 1){
                                print join("\t", ($chr, $i, "-", $type, $meth{$chr}->{"-"}->{$i}, $cover{$chr . "-"}->[$i])), "\n";
                                if($meth{$chr}->{"-"}->{$i} >= $meth_percent_cutoff * $cover{$chr . "-"}->[$i]){
                                        $num_meth_sites{$chr . "-"}->{$type}++;
                                }
                        }
                }
        $pos++;
        }
}
if($debug){
        foreach my $type("CG", "CHG", "CHH"){
        foreach my $chr(sort keys %chr_len){
          print STDERR "$chr\t$type\t", $num_meth_sites{$chr . "-"}->{$type}, "\n";
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
                if(!defined $num_meth_sites{$chr . $str}->{$type}){
                        $num_meth_sites{$chr . $str}->{$type} = 0;
                }
                $total_C += $num_sites{$chr . $str}->{$type};
                $total_mC += $num_meth_sites{$chr . $str}->{$type};
                print STDERR join("\t", ($str, $chr, $type, $num_sites{$chr . $str}->{$type}, $num_meth_sites{$chr . $str}->{$type}
                ,$num_meth_sites{$chr . $str}->{$type}/$num_sites{$chr . $str}->{$type}
                )), "\n";
        }
        print STDERR join("\t", ($str, "AllChrs", $type, $total_C, $total_mC, $total_mC/$total_C)), "\n"; 
}
}

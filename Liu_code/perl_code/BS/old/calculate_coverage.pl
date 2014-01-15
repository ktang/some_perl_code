#!/usr/bin/perl -w
# given BSMAP output, calculate the depth of coverage at each genome position
# and percentage of genome that was covered by at least one or two reads.
use strict;
my %cover;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478);
while(<>){
        chomp;
        my ($id, $seq, $map_flag, $chr, $start, $strand, $num_mm, $mm_info, $mC_loc)        = split /\t/;
        foreach my $i($start..($start+length($seq)-1)){
                #$cover{$chr}->[$i] = $cover{$chr}->[$i] + 1;
                if($strand eq "++" || $strand eq "+-"){
                    $cover{$chr . "W"}->[$i]++;
                }else{
                        $cover{$chr . "C"}->[$i]++;
                }
        }
}
my (%cover_1W, %cover_2W, %cover_1C, %cover_2C);
my $total_len = 0;
foreach my $chr(keys %chr_len){
        $total_len += $chr_len{$chr};
        foreach my $i(1..$chr_len{$chr}){
                if(defined $cover{$chr . "W"}->[$i]){ 
                        if($cover{$chr . "W"}->[$i] >= 1){
                            $cover_1W{$chr}++;
                            $cover_1W{"total"}++;
                        }
                        if($cover{$chr . "W"}->[$i] >= 2){
                                $cover_2W{$chr}++;
                                $cover_2W{"total"}++;
                        }
                }
                if(defined $cover{$chr . "C"}->[$i]){
                        if($cover{$chr . "C"}->[$i] >= 1){
                            $cover_1C{$chr}++;
                            $cover_1C{"total"}++;
                        }
                        if($cover{$chr . "C"}->[$i] >= 2){
                                $cover_2C{$chr}++;
                                $cover_2C{"total"}++;
                        }
                }
        }
}
print "Coverage >= 1:\n";
print "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
        print $chr, "\t", $cover_1W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_1W{"total"}/$total_len * 100, "%\n";
print "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
        print $chr, "\t", $cover_1C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_1C{"total"}/$total_len * 100, "%\n";
print "Coverage >= 2:\n";
print "On Watson strand:\n";
foreach my $chr(sort keys %chr_len){
        print $chr, "\t", $cover_2W{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_2W{"total"}/$total_len * 100, "%\n";
print "On Crick strand:\n";
foreach my $chr(sort keys %chr_len){
        print $chr, "\t", $cover_2C{$chr}/$chr_len{$chr} * 100, "%\n";
}
print "Whole genome", "\t", $cover_2C{"total"}/$total_len * 100, "%\n";  

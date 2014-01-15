#!/usr/bin/perl -w
# parse BS_Seeker outputs and produce a bed file that indicates each methylated sites
use strict;
my $print_counts = 1;
my %chrs = ('0001'=>'chr1','0002'=>'chr2','0003'=>'chr3','0004'=>'chr4',
            '0005'=>'chr5', '0006'=>'chrM', '0007'=>'chrC');
my %meth; # chr->location= [type, #reads methylated, # reads unmethylated];
while(<>){
        chomp;
        my ($read, $mism, $strand, $loc, $refSeq, $readSeq, $methyCode) = split /\t/;
        my ($chrNum, $startLoc) = (0,0);
        if($loc =~ /(\d+)[\+|\-](\d+)/){
                ($chrNum, $startLoc) = ($1, $2);
        }else{
                die "read $read match location $loc";
        }
        my $chr = $chrs{$chrNum};
        my @chars = split ("", $methyCode);
        foreach my $i(0..$#chars){
                if($chars[$i] ne '-'){
                        my $thisLoc = $startLoc+$i;
                        if(defined $meth{$chr}->{$thisLoc}){
                                if($meth{$chr}->{$thisLoc}->[0] eq (uc $chars[$i])){
                                        if($meth{$chr}->{$thisLoc}->[0] eq $chars[$i]){
                                                $meth{$chr}->{$thisLoc}->[1]++;
                                        }else{
                                                $meth{$chr}->{$thisLoc}->[2]++;
                                        }
                                }else{
                                        die "read $read, $loc, has different type of methylation";
                                }
                        }else{
                                if($chars[$i] eq 'X' || $chars[$i] eq 'Y' || $chars[$i] eq 'Z'){
                                        $meth{$chr}->{$thisLoc}->[0] = $chars[$i];
                                        $meth{$chr}->{$thisLoc}->[1] = 1;
                                }else{
                                        $meth{$chr}->{$thisLoc}->[0] = uc $chars[$i];
                                        $meth{$chr}->{$thisLoc}->[2] = 1;
                                }
                        }
                }
        }
}

my %scores = ('X' => 2, 'Y'=>1.5, 'Z'=>1);
my %names = ('X'=>'CG', 'Y'=>'CHG', 'Z'=>'CHH');
foreach my $chr(sort keys %meth){
        my %locs = %{$meth{$chr}};
        foreach my $loc(sort {$a <=>$b} keys %locs){
                print join("\t", ($chr, $loc, $loc, $names{$locs{$loc}->[0]}, $scores{$locs{$loc}->[0]}));
                if($print_counts){
                        print "\t", $locs{$loc}->[1], "\t", $locs{$loc}->[2];
                }
                print "\n";
        }
}


#!/usr/bin/perl -w
# #compare two wig files and identify hyper/hypo methylated regions
#Kai input is .gff
#/Users/tang58/Nimblegen_project/Ecker_data/mc_col0_tair9_no_warn.gff
#/Users/tang58/Nimblegen_project/Ecker_data/mc_rdd_tair9_no_warn.gff
use strict;

my $usage = "$0 <win_size> <step_size> <min_fold_change>";
die $usage unless (@ARGV ==3 );
#my ($wig1, $wig2, $isHypo) = @ARGV[0..2];

#my $window_size = 200;
#my $step_size = 100;
#my $min_fold_change = 3;
my ($window_size, $step_size, $min_fold_change) = @ARGV[0..2];

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my $rdd_file = "/Users/tang58/Nimblegen_project/Ecker_data/mc_rdd_tair9_no_warn.gff";
my $col_file = "/Users/tang58/Nimblegen_project/Ecker_data/mc_col0_tair9_no_warn.gff";

#my (%meth1, %meth2);
my (%density1, %density2);# 1 for col0 and 2 for rdd
#open(MF1, $wig1) or die "Can't open $wig1:$!";
open(COL, $col_file) or die "Can't open $col_file:$!";
open(RDD, $rdd_file) or die "Cannot open $rdd_file:$!";

my $chr = "";
my $num_windows = int($window_size / $step_size + 0.99);

#while(<MF1>){
while(<COL>){
#	next if(/^track/);
#	if(/variableStep\s+chrom=(\w+)/){
#	    $chr = $1;
#		next;
#	}
	chomp;
	#my ($pos, $val) = split /\t/;
	#$meth1{$chr}->{$pos} = $val;
	my @a = split "\t";
	my $pos = $a[3];
	my $chr = lc $a[0];
	my $base_index = int($pos / $step_size);
	foreach my $i(($base_index-$num_windows)..($base_index+$num_windows)){#for each postion, see all the potentional window,to
	#check wether it in the window. smarter than my algrothim.
		if($i >= 0){
			if($pos >= 1 && $pos >= ($i*$step_size) && $pos < ($i*$step_size+$window_size)){
				$density1{$chr}->[$i]++;
			}
		}
	}
}

#open(MF2, $wig2) or die "Can't open $wig2:$!";
#while(<MF2>){
#	next if(/^track/);
#	if(/variableStep\s+chrom=(\w+)/){
#	    $chr = $1;
#		next;
#	}
while(<RDD>){
	chomp;
	#my ($pos, $val) = split /\t/;
	#$meth2{$chr}->{$pos} = $val;
	my @a = split "\t";
	my $chr = lc $a[0];
	my $pos = $a[3]; 
    my $base_index = int($pos / $step_size);
	foreach my $i(($base_index-$num_windows)..($base_index+$num_windows)){
		if($i >= 0){
			if($pos >= 1 && $pos >= ($i*$step_size) && $pos < ($i*$step_size+$window_size)){
				$density2{$chr}->[$i]++;
			}
		}
	}

}

#print join("\t", ("Chr", "window_center", $ARGV[0], $ARGV[1], "difference")), "\n";
#my $desc = "hypo-methylated";
#my $color = "0,250,0";
#if($isHypo == 0){
#	$desc = "hyper-methylated";
#	$color = "250,0,0";
#}
#print "track name=", $desc, 'description="', $desc, ' regions" color=', $color, ' useScore=1', "\n";

my ($accum1, $accum2, $last_start, $last_end, $last_chr) = (0, 0, 0, 0, "NONE");
foreach my $chr(sort keys %chr_len){
	my $len = $chr_len{$chr};
	foreach my $i(0..(int($len/$step_size)+1)){
		if(!defined $density1{$chr}->[$i]){
			$density1{$chr}->[$i] = 0;
		}
		if(!defined $density2{$chr}->[$i]){
			$density2{$chr}->[$i] = 0;
		}
		
=head		
		my $hypo_change = $density1{$chr}->[$i];
		if($density2{$chr}->[$i] > 0){
			$hypo_change = $density1{$chr}->[$i]/$density2{$chr}->[$i];
		}
		if($hypo_change >= $min_fold_change){
			if($isHypo && $density1{$chr}->[$i] >= $min_fold_change * 2){
				my ($new_start, $new_end) = ($i*$step_size+1, $i*$step_size+$window_size);
				if($chr ne $last_chr || $new_start > $last_end + 5){ # not overlap
					my $fold = $accum1;
					if($accum2 > 0){
						$fold = $accum1/$accum2;
					}
					if($fold > 10){
						$fold = 10;
					}
					my $score = int($fold*100);
					
					print join("\t", ($last_chr, $last_start, $last_end, ".",$score)), "\n";
					($last_chr, $last_start, $last_end, $accum1, $accum2) = (
					 $chr, $new_start, $new_end, $density1{$chr}->[$i], $density2{$chr}->[$i]);
				}else{
					$last_end = $i*$step_size+$window_size;
					$accum1 += $density1{$chr}->[$i];
					$accum2 += $density2{$chr}->[$i];
				}
			}
		}
=cut		
		my $hyper_change = $density2{$chr}->[$i];
		if($density1{$chr}->[$i] > 0){
			$hyper_change = $density2{$chr}->[$i]/$density1{$chr}->[$i];
		}

		if($hyper_change >= $min_fold_change){
        #   if(!$isHypo && $density2{$chr}->[$i] >= $min_fold_change * 2){
			if($density2{$chr}->[$i] >= $min_fold_change * 2){
	            my ($new_start, $new_end) = ($i*$step_size+1, $i*$step_size+$window_size);
				if($chr ne $last_chr || $new_start > $last_end + 5){ # not overlap
					my $fold = $accum2;
					if($accum1 > 0){
						$fold = $accum2/$accum1;
					}
					if($fold > 10){
						$fold = 10;
					}
					my $score = int($fold*100);
					
					print join("\t", ($last_chr, $last_start, $last_end, ".",$score)), "\n";
					($last_chr, $last_start, $last_end, $accum1, $accum2) = (
					 $chr, $new_start, $new_end, $density1{$chr}->[$i], $density2{$chr}->[$i]);
				}else{
					$last_end = $i*$step_size+$window_size;
					$accum1 += $density1{$chr}->[$i];
					$accum2 += $density2{$chr}->[$i];
				}
			}
		}
	}
}
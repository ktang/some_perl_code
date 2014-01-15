#!/usr/bin/perl -w
# two soap output files get enrich regions

use strict;

my $usage = "$0 <ctrl_soap_file> <IP_file> <win_size> <step_size> <min_fold_change> <min_IP_reads> <output>";
die $usage unless (@ARGV == 7);

my ($ctrl_file, $IP_file, $window_size, $step_size, $min_fold_change) = @ARGV[0..4]; 


my $min_IP_reads = $ARGV[5]; #10;

my $output = $ARGV[6];

die "$output exists!!" unless !(-e $output) ; 

open (OUT,">$output") or die "cannot open $output";

print OUT join ("\t", ("Chr", "Start", "End", "accum_IP", "accum_ctrl", "score")), "\n";
#my $window_size = 200;
#my $step_size = 100;
#my $min_fold_change = 3;

my %chr_len = ("Chr1"=>30427671, "Chr2"=>19698289, "Chr3"=>23459830, "Chr4"=>18585056, "Chr5"=>26975502);

#my ($wig1, $wig2, $isHypo) = @ARGV[0..2];
#my (%meth1, %meth2);
#my (%density1, %density2);

my (%den_ctrl, %den_IP);

#open(MF1, $wig1) or die "Can't open $wig1:$!";

open (CTRL, $ctrl_file) or die "cannot open $ctrl_file";
open (IP, $IP_file) or die "cannot open $IP_file";

my $chr = "";

my $num_windows = int($window_size / $step_size + 0.99);
#while(<MF1>){
#	next if(/^track/);
#	if(/variableStep\s+chrom=(\w+)/){
#	    $chr = $1;
#		next;
#	}

#	my ($pos, $val) = split /\t/;
	#$meth1{$chr}->{$pos} = $val;
while(<CTRL>){
	chomp;
	my @a = split "\t";
	my $pos = $a[8];
	my $chr = $a[7];
	my $base_index = int($pos / $step_size);
	foreach my $i(($base_index-$num_windows)..($base_index+$num_windows)){
		if($i >= 0){
			if($pos >= ($i*$step_size) && $pos < ($i*$step_size+$window_size)){
				$den_ctrl{$chr}->[$i]++;
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
while (<IP>){
	chomp;
#	my ($pos, $val) = split /\t/;
	#$meth2{$chr}->{$pos} = $val;
	my @a = split "\t";
	my $pos = $a[8];
	my $chr = $a[7];

    my $base_index = int($pos / $step_size);
	foreach my $i(($base_index-$num_windows)..($base_index+$num_windows)){
		if($i >= 0){
			if($pos >= ($i*$step_size) && $pos < ($i*$step_size+$window_size)){
				$den_IP{$chr}->[$i]++;
			}
		}
	}

}

#print join("\t", ("Chr", "window_center", $ARGV[0], $ARGV[1], "difference")), "\n";

my ($accum_IP, $accum_ctrl, $last_start, $last_end, $last_chr) = (0, 0, 0, 0, "NONE");

foreach my $chr(sort keys %chr_len){
	my $len = $chr_len{$chr};
	foreach my $i(0..(int($len/$step_size)+1)){
		if(!defined $den_IP{$chr}->[$i]){
			$den_IP{$chr}->[$i] = 0;
		}
		if(!defined $den_ctrl{$chr}->[$i]){
			$den_ctrl{$chr}->[$i] = 0;
		}
		
		my $change = $den_IP{$chr}->[$i];
		
		if($den_ctrl{$chr}->[$i] > 0){
			$change = $den_IP{$chr}->[$i]/$den_ctrl{$chr}->[$i];
		}
		
		if($change >= $min_fold_change){
			if($den_IP{$chr}->[$i] >= $min_IP_reads){#$min_fold_change * 2)
				
				my ($new_start, $new_end) = ($i*$step_size+1, $i*$step_size+$window_size);
				
				if($chr ne $last_chr || $new_start > $last_end + 5){ # not overlap
					my $fold = $accum_IP;
					if($accum_ctrl > 0){
						$fold = $accum_IP/$accum_ctrl;
					}
					if($fold > 10){
						$fold = 10;
					}
					my $score = int($fold*100);
					
					print OUT join("\t", ($last_chr, $last_start, $last_end, $accum_IP,$accum_ctrl,$score)), "\n";
					($last_chr, $last_start, $last_end, $accum_IP, $accum_ctrl) = (
						 $chr, $new_start, $new_end, $den_IP{$chr}->[$i], $den_ctrl{$chr}->[$i]);
				}else{
					$last_end = $i*$step_size+$window_size;
					$accum_IP += $den_IP{$chr}->[$i];
					$accum_ctrl += $den_ctrl{$chr}->[$i];
				}
			}
		}
	}
}


#	my $hyper_change = $den_ctrl{$chr}->[$i];
#		if($den_IP{$chr}->[$i] > 0){
#			$hyper_change = $den_ctrl{$chr}->[$i]/$den_IP{$chr}->[$i];
#		}

#		if($hyper_change >= $min_fold_change){
#           if(!$isHypo && $den_ctrl{$chr}->[$i] >= $min_fold_change * 2){
#	            my ($new_start, $new_end) = ($i*$step_size+1, $i*$step_size+$window_size);
#				if($chr ne $last_chr || $new_start > $last_end + 5){ # not overlap
#					my $fold = $accum_ctrl;
#					if($accum_IP > 0){
#						$fold = $accum_ctrl/$accum_IP;
#					}
#					if($fold > 10){
#						$fold = 10;
#					}
#					my $score = int($fold*100);
#					
#					print join("\t", ($last_chr, $last_start, $last_end, ".",$score)), "\n";
#					($last_chr, $last_start, $last_end, $accum_IP, $accum_ctrl) = (
#					 $chr, $new_start, $new_end, $den_IP{$chr}->[$i], $den_ctrl{$chr}->[$i]);
#			}else{
#				$last_end = $i*$step_size+$window_size;
#				$accum_IP += $den_IP{$chr}->[$i];
#				$accum_ctrl += $den_ctrl{$chr}->[$i];
#			}
#		 }
#	  }


#	}
#}
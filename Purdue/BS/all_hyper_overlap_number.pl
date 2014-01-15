#!/usr/bin/perl -w
#strand specific 

use strict;

my $debug = 0;

#my $window_size = 200;
#my $step_size = 100;
my $min_fold = 2;
my $min_cutoff = 10; #mC_ros3 at least >= $min_cutoff

#my $wig_Col0 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Col0_s_1_methylation.wig";
#my $wig_ros3 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Ros3_2_s_8_methylation.wig";
#my $wig_Col0 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_col0_C.sgr";
#my $wig_ros3 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_rdd_C.sgr";

#if($debug){
#	$wig_Col0 = "/Users/tang58/try/Ros3_binding/wig_files/Col0.wig";
#	$wig_ros3 = "/Users/tang58/try/Ros3_binding/wig_files/Ros3.wig";	
#}

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

#my $usage = "$0 <indir> <outdir> <rdd|ros1|ros3> <fold> <min_number>";


#die $usage unless (@ARGV == 5);

#my ($indir, $outdir, $label) = @ARGV[0..2];

#my ($min_fold_change, $min_cutoff) = @ARGV[3..4];

#die "wrong directory\n\n" unless (-e $indir and -e $outdir);

#die "rdd|ros1|ros3\n\n" unless ($label eq "rdd" or $label eq "ros1" or $label eq "ros3"
#	or $label eq "hsp" or $label eq "phd" or $label eq "zdp" or $label eq "ros2" or $label eq "ros4" );

my $wig_col0 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ1_Col0.wig";
my $ros2_wig = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ6_ros2.wig";

my @files = ("/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_17_rdd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_18_rdd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ3_ros1_4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ4_ros1_4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ7_ros3_2.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ8_ros3_2.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ9_ros4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ10_ros4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_15_zdp.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_16_zdp.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ12_phd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ13_phd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ11_hsp.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ12_hsp.wig"
			 );

if ($debug) {print STDERR join("\n", @files), "\n"}

my (%col0_h, %ros2_h);

read_first_file($wig_col0, \%col0_h);
read_first_file($ros2_wig, \%ros2_h);

my (%ros1_h, %ros3_h, %ros4_h, %rdd_h, %zdp_h, %phd_h, %hsp_h);
my (%ros1_temp, %ros3_temp, %ros4_temp, %rdd_temp, %zdp_temp, %phd_temp, %hsp_temp);

read_first_file($files[0], \%rdd_temp);
read_first_file($files[2], \%ros1_temp);
read_first_file($files[4], \%ros3_temp);
read_first_file($files[6], \%ros4_temp);
read_first_file($files[8], \%zdp_temp);
read_first_file($files[10], \%phd_temp);
read_first_file($files[12], \%hsp_temp);

read_second_file($files[1], \%rdd_h, \%rdd_temp);
read_second_file($files[3], \%ros1_h, \%ros1_temp);
read_second_file($files[5], \%ros3_h, \%ros3_temp);
read_second_file($files[7], \%ros4_h, \%ros4_temp);
read_second_file($files[9], \%zdp_h, \%zdp_temp);
read_second_file($files[11], \%phd_h, \%phd_temp);
read_second_file($files[13], \%hsp_h, \%hsp_temp);

#0:rdd, 1-4: ros1-4, 5: zdp,  6:phd,  7:hsp
my @nums;
for my $i (0..7){
	for my $j(0..7){
		$nums[$i][$j] = 0;
	}
}

my $outdir = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/1kb_hyper";

if (-e "$outdir/rdd_hyper.bed") {die "output exists\n"}
open (RDD, ">$outdir/rdd_hyper.bed"); 
open (ROS1, ">$outdir/ros1_hyper.bed"); 
open (ROS2, ">$outdir/ros2_hyper.bed"); 
open (ROS3, ">$outdir/ros3_hyper.bed"); 
open (ROS4, ">$outdir/ros4_hyper.bed"); 
open (ZDP, ">$outdir/zdp_hyper.bed"); 
open (PHD, ">$outdir/phd_hyper.bed"); 
open (HSP, ">$outdir/hsp_hyper.bed"); 

foreach my $chr (sort keys %chr_len){
	my $total_index = int($chr_len{$chr} / 1000) + 1;
	for my $i (0..$total_index - 1){
		my $num_col0 = count($chr, $i, \%col0_h);
		my $num_rdd = count($chr, $i, \%rdd_h);
		my $num_ros1 = count($chr, $i, \%ros1_h);
		my $num_ros2 = count($chr, $i, \%ros2_h);
		my $num_ros3 = count($chr, $i, \%ros3_h);
		my $num_ros4 = count($chr, $i, \%ros4_h);
		my $num_zdp = count($chr, $i, \%zdp_h);
		my $num_phd = count($chr, $i, \%phd_h);
		my $num_hsp = count($chr, $i, \%hsp_h);
		my @flags = (0, 0, 0, 0, 0, 0, 0, 0);
		if ($num_col0 == 0){
			if ($num_rdd >= $min_cutoff){print RDD join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_rdd , $num_col0 ) ),"\n" ; $nums[0][0]++; $flags[0] = 1;  }
			if ($num_ros1 >= $min_cutoff){print ROS1 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros1 , $num_col0  ) ),"\n" ; $nums[ 1][ 1]++; $flags[1] = 1  }
			if ($num_ros2 >= $min_cutoff){print ROS2 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros2 , $num_col0  ) ),"\n" ; $nums[ 2][2]++; $flags[2] = 1 }
			if ($num_ros3 >= $min_cutoff){print ROS3 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros3 , $num_col0  ) ),"\n" ; $nums[ 3][3]++; $flags[3] = 1 }
			if ($num_ros4 >= $min_cutoff){print ROS4 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros4 , $num_col0  ) ),"\n" ; $nums[ 4][4]++; $flags[4] = 1 }
			if ($num_zdp >= $min_cutoff){print ZDP join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_zdp , $num_col0  ) ),"\n" ; $nums[ 5][ 5]++; $flags[5] = 1 }
			if ($num_phd >= $min_cutoff){print PHD join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_phd , $num_col0  ) ),"\n" ; $nums[ 6][6]++; $flags[6] = 1 }
			if ($num_hsp >= $min_cutoff){print HSP join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_hsp , $num_col0  ) ),"\n" ; $nums[ 7][7]++; $flags[7] = 1 }
		}else{
			if ($num_rdd >= $min_cutoff and $num_rdd/$num_col0 >= $min_fold ){print RDD join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_rdd , $num_col0  ) ),"\n" ; $nums[0][ 0]++; $flags[0] = 1 }
			if ($num_ros1 >= $min_cutoff and $num_ros1/$num_col0 >= $min_fold ){print ROS1 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros1 , $num_col0  ) ),"\n" ; $nums[ 1][ 1]++; $flags[1] = 1  }
			if ($num_ros2 >= $min_cutoff and $num_ros2/$num_col0 >= $min_fold ){print ROS2 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros2 , $num_col0  ) ),"\n" ; $nums[ 2][2]++; $flags[2] = 1 }
			if ($num_ros3 >= $min_cutoff and $num_ros3/$num_col0 >= $min_fold ){print ROS3 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros3 , $num_col0  ) ),"\n" ; $nums[ 3][3]++; $flags[3] = 1 }
			if ($num_ros4 >= $min_cutoff and $num_ros4/$num_col0 >= $min_fold ){print ROS4 join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_ros4 , $num_col0  ) ),"\n" ; $nums[ 4][4]++; $flags[4] = 1 }
			if ($num_zdp >= $min_cutoff and $num_zdp/$num_col0 >= $min_fold ){print ZDP join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_zdp , $num_col0  ) ),"\n" ; $nums[ 5][ 5]++; $flags[5] = 1 }
			if ($num_phd >= $min_cutoff and $num_phd/$num_col0 >= $min_fold ){print PHD join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_phd , $num_col0  ) ),"\n" ; $nums[ 6][ 6]++; $flags[6] = 1 }
			if ($num_hsp >= $min_cutoff and $num_hsp/$num_col0 >= $min_fold ){print HSP join ("\t", ( $chr, 1000 * $i + 1, 1000 * ($i + 1), $num_hsp , $num_col0  ) ),"\n" ; $nums[ 7][7]++; $flags[7] = 1 }
		}
		
		for my $i(0..6){
			if ($flags[$i] == 1){
				for my $j ($i + 1 .. 7){
					if($flags[$j] == 1){
						$nums[$i][$j] ++;
						$nums[$j][$i] ++;

					}
				}
			}
		}		
	}
}

open (OUT, ">$outdir/hyper_preliminary_result.txt");

print OUT join ("\t", ( "sample", "rdd", "ros1", "ros2", "ros3", "ros4", "zdp", "phd", "hsp")), "\n";
print OUT join ("\t", ("rdd", @{$nums[0]}[0..7] ) ), "\n";
print OUT join ("\t", ("ros1", @{$nums[1]}[0..7] ) ), "\n";
print OUT join ("\t", ("ros2", @{$nums[2]}[0..7]) ), "\n";
print OUT join ("\t", ("ros3", @{$nums[3]}[0..7]) ), "\n";
print OUT join ("\t", ("ros4", @{$nums[4]}[0..7]) ), "\n";
print OUT join ("\t", ("zdp", @{$nums[5]}[0..7]) ), "\n";
print OUT join ("\t", ("phd", @{$nums[6]}[0..7]) ), "\n";
print OUT join ("\t", ("hsp", @{$nums[7]}[0..7]) ), "\n";

exit;

sub read_first_file{
	my ($input, $hash_ref) = @_;
	print STDERR "reading $input...\n";
	open (IN, $input) or die "cannot open $input:$!";
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
		${$hash_ref}{$chr}->[$pos] = $val;
	}	
	print STDERR "$input done\n";
}

sub read_second_file{
	my ($input, $hash_ref, $former) = @_;
	print STDERR "reading $input...\n";
	open (IN, $input) or die "cannot open $input:$!";
	my $chr = "";
	while (<IN>){
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = lc $1;
			next;
		}
		chomp;
		my ($pos, $val) = split /\t/;
		if (defined ${$former}{$chr}->[$pos]){
			${$hash_ref}{$chr}->[$pos] = $val;
		}
	}
	print STDERR "$input done\n";
}

sub count {
	my ($chr, $index, $hashref) = @_;
	my $number = 0;
	for my $i (1000 * $index + 1 .. 1000 * ($index + 1)){
		if (defined ${$hashref}{$chr}->[$i]){
			$number++;
		}
	}
	return $number;
}
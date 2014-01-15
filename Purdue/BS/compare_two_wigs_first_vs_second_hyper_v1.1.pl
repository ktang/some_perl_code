#!/usr/bin/perl -w
# use moving window compare two mutants(including WT).
# the last parameter is and/or, which indicate how to 
#  use biological replicate. 

# v1.1 parameter input from cmd

use strict;

my $debug = 0;

#my $window_size = 200;
#my $step_size = 100;
#my $min_fold_change = 3;
#my $min_cutoff = 10; #mC_ros3 at least >= $min_cutoff

#my $wig_Col0 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Col0_s_1_methylation.wig";
#my $wig_ros3 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Ros3_2_s_8_methylation.wig";
#my $wig_Col0 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_col0_C.sgr";
#my $wig_ros3 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_rdd_C.sgr";

#if($debug){
#	$wig_Col0 = "/Users/tang58/try/Ros3_binding/wig_files/Col0.wig";
#	$wig_ros3 = "/Users/tang58/try/Ros3_binding/wig_files/Ros3.wig";	
#}

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);
my %labels = (ros1=>1, ros2=>1, ros3=>1, ros4=>1, rdd=>1, zdp=>1, phd=>1, hsp=>1, col0=>1);

my $usage = "$0 <label1> <label2> <outdir> <output_pre> <and/or> <WinSize> <step_size> <min_fold_change> <min_cutoff>";

die $usage unless (@ARGV == 9);

my ($window_size, $step_size, $min_fold_change, $min_cutoff ) = @ARGV[5..8];



my ($label1, $label2, $outdir, $output_pre) = @ARGV[0..3];

my $logic = $ARGV[4];

die unless ($logic eq "and" or $logic eq "or");

die "wrong directory\n\n" unless (-e $outdir);
die "wrong label\n" unless (defined $labels{$label1} and defined $labels{$label2} and $label1 ne $label2);

my $output = $output_pre."_$logic"."_Size$window_size"."Step$step_size"."Fold$min_fold_change"."Cutoff$min_cutoff".".txt";

print STDERR "output:$output\n";

if (-e "$outdir/$output") {die "output exists!!"}


my $num_windows = int($window_size / $step_size + 0.99);


open (OUT, ">$outdir/$output") or die "cannot open $output:$!";

print OUT join("\t", ("chr", "start", "end", $label1."_score",  $label1."_num",
			 $label2."_score",  $label2."_num")), "\n";

my ($file_first1, $file_first2, $file_second1, $file_second2);

#my $wig_col0 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ1_Col0.wig";
my $wig_col0 = "/Volumes/Macintosh_HD_2/Kai_BS/brat_output/wig/col0_s_8_JKZ131.wig";

my @files = ("/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_17_rdd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_18_rdd.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ3_ros1_4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ4_ros1_4.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ5_ros2.wig",
			 "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ6_ros2.wig",
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


if ($label1 eq "col0"){
	$file_first1 = $wig_col0;
	$file_first2 = $wig_col0;
}
else{
	($file_first1, $file_first2) = grep /$label1/, @files;
}

if ($label2 eq "col0"){
	$file_second1 = $wig_col0;
	$file_second2 = $wig_col0;
}
else{
	($file_second1, $file_second2) = grep /$label2/, @files;
}

print STDERR join("\n",($file_first1, $file_first2, $file_second1, $file_second2)),"\n";



my (%first_hh, %second_hh);#hh means hash_hash
my (%first0_hh, %second0_hh);#hh means hash_hash

if ( $logic eq "and"){
	
	read_first_file_hash($file_first1, \%first0_hh);
	print STDERR "$label1 met number:and_first:", cal_met(\%first0_hh) , "\n";
	read_second_file_hash($file_first2, \%first_hh, \%first0_hh);
	print STDERR "$label1 met number:and_second:", cal_met(\%first_hh) , "\n";
		
	read_first_file_hash($file_second1, \%second0_hh);
	print STDERR "$label2 met number:and_second:", cal_met(\%second0_hh) , "\n";
	read_second_file_hash($file_second2, \%second_hh, \%second0_hh);
	print STDERR "$label2 met number:and_second:", cal_met(\%second_hh) , "\n";
}
elsif($logic eq "or"){
	
	read_first_file_hash($file_first1, \%first_hh);
	print STDERR "$label1 met number:or_first:", cal_met(\%first_hh) , "\n";
	read_first_file_hash($file_first2, \%first_hh);
	print STDERR "$label1 met number:or_first:", cal_met(\%first_hh) , "\n";
	
	read_first_file_hash($file_second1, \%second_hh);
	print STDERR "$label2 met number:or_second:", cal_met(\%second_hh) , "\n";
	read_first_file_hash($file_second2, \%second_hh);	
	print STDERR "$label2 met number:or_second:", cal_met(\%second_hh) , "\n";
}
else{
	die "wrong logic\n";
}


#print STDERR "$label1 met number:and_first", cal_met(\%first0_hh) , "\n";
#print STDERR "$label2 met number:and_first", cal_met(\%first0_hh) , "\n";

#print STDERR "$label1 met number:", cal_met(\%first_hh) , "\n";
#print STDERR "$label2 met number:", cal_met(\%second_hh) , "\n";


my (%density1, %density2);

cal_density(\%first_hh, \%density1);
cal_density(\%second_hh, \%density2);



my ($last_start, $last_end)  = (0, 0);
my $last_chr = "NONE";


foreach my $chr(sort keys %chr_len){
	my $len = $chr_len{$chr};
	
	foreach my $i(0..(int($len/$step_size)+1)){
		if(!defined $density1{$chr}->[$i]){
			$density1{$chr}->[$i] = 0;
		}
		if(!defined $density2{$chr}->[$i]){
			$density2{$chr}->[$i] = 0;
		}
		
		my $change = $density1{$chr}->[$i];
		if($density2{$chr}->[$i] > 0){
			$change = $density1{$chr}->[$i] / $density2{$chr}->[$i];
		}
		
		if ( $change >=  $min_fold_change and $density1{$chr}->[$i] >= $min_cutoff){
			my ($new_start, $new_end) = ($i*$step_size+1, $i*$step_size+$window_size);
			
			if($chr ne $last_chr || $new_start > $last_end + 5){ # not overlap
				
				my ($score1, $num1, $score2, $num2) = (0, 0, 0, 0);
			
				foreach my $k ($last_start..$last_end){
					if(defined $first_hh{$last_chr}->{$k}){$num1++; $score1 += abs($first_hh{$last_chr}->{$k}) }
					if(defined $second_hh{$last_chr}->{$k}){$num2++; $score2 += abs($second_hh{$last_chr}->{$k}) }
				}
			
				print OUT join("\t", ($last_chr, $last_start, $last_end, $score1, $num1, $score2, $num2)), "\n" unless ($last_chr eq "NONE");
				($last_chr, $last_start, $last_end) = ($chr, $new_start, $new_end);
			}
			else{
				$last_end = $i*$step_size+$window_size;
			}
		}
	}
}

close(OUT);
exit;

sub cal_density{
	my ($ref_file, $ref_den) = @_;
	
	print STDERR "cal_density...\t";
	
	
	foreach my $chr (sort keys %{$ref_file} ){
	
		foreach my $pos (sort {$a <=> $b} keys %{ ${$ref_file}{$chr}} ){
			
			my $base_index = int($pos / $step_size);
			
			foreach my $i(($base_index-$num_windows)..($base_index+$num_windows)){
				if($i >= 0){
					if($pos >= 1 and $pos > ($i*$step_size) and $pos <= ($i*$step_size+$window_size)){
						${$ref_den}{$chr}->[$i]++;
					}
				}
			}
		}
	}
	print STDERR "DONE\n";
}

sub read_first_file_hash{
	my ($input, $hash_ref) = @_;
	print STDERR "read_first_file_hash reading first $input...\t";
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
		if(defined ${$hash_ref}{$chr}->{$pos}){
			my $temp = (${$hash_ref}{$chr}->{$pos} + $val) / 2.0;
			${$hash_ref}{$chr}->{$pos} = $temp;
		}
		else{
			${$hash_ref}{$chr}->{$pos} = $val;
		}	
	}	
	close(IN);
	print STDERR "DONE\n";
}

sub read_second_file_hash{
	my ($input, $hash_ref, $former) = @_;
	print STDERR "read_second_file_hash reading $input...\t";
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
		if (defined ${$former}{$chr}->{$pos}){
			my $temp = (${$former}{$chr}->{$pos} + $val) / 2.0;
			${$hash_ref}{$chr}->{$pos} = $temp;
			
		}
	}
	close(IN);
	print STDERR "DONE\n";
}

sub cal_met{
	my ($ref) = @_;
	my $num = 0;
	foreach my $chr (sort keys %{$ref}){
		$num += scalar(keys %{${$ref}{$chr}})
	}
	return $num;
}
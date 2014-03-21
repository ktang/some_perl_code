#!/usr/bin/perl -w
# simulation: randomly generate the same length of the ros1-4 target loci,

#Chr output big letter used for histone sum

use strict;

my $debug = 0;
if($debug){
	print STDERR "debug: $debug\n\n";
}

my $usage = "$0 \n <input_bed_like_DMR> <output> \n\n";
die $usage unless(@ARGV == 2);

#my $itr_times = shift or die "itr";
#my $wig_file = shift or die "wig";
my  $target_file = shift or die ;
my $output = shift or die "output";

#die unless (-e $wig_file);
#print STDERR "wig file used:\n$wig_file\n\n";

die if (-e $output);

#open(IN, $input) or die "cannot open $input: $!";
open(OUT, ">$output") or die "cannot open $output: $!";

#my $target_file = "/Volumes/Macintosh_HD_2/idm1_new_met_data/script/Liu_algorithm_Oct1_data_set/output_both_depth100_Oct3/1000_1000_5_10_P001/raw_with_head/ros1_4_hyper_P0.01_reduced_boundary_both_depth100_WinSize1000_gap1000_initialCutoff5_reportCutoff10.txt";
die "target_file" unless (-e $target_file);

my @loci_lengths;
my $total_length;
read_loci_list($target_file, \@loci_lengths, \$total_length);


my @regions = ();
	
generate_loci(\@regions);

print OUT join("\t", ("chr", "start", "end")), "\n";	
for my $k (0..$#loci_lengths){
	my ($chr, $start, $end) = @{$regions[$k]};
	$chr =~ s/c/C/;
	print OUT join("\t", ( $chr, $start, $end )), "\n";
	
}
close(OUT);
exit;

sub generate_loci{
	my ($ref) = @_;
	my ($len_chr1,$len_chr2,$len_chr3,$len_chr4,$len_chr5)=
		(30427671, 19698289, 23459830, 18585056, 26975502);
#	my $len_chr1  = 30427671;
	my $point_2   = 50125960;   # 30427671 + 19698289,
	my $point_3   = 73585790;   # 30427671 + 19698289 + 23459830
	my $point_4   = 92170846;   # 30427671 + 19698289 + 23459830 +  18585056
	my $total_chr = 119146348;  # 30427671 + 19698289 + 23459830 +  18585056 + 26975502
	
#	if($debug){
#		print STDERR "in generate_loci: loci_lengths = ", $#loci_lengths, "\n\n";
#	}
	
	foreach my $i (0..$#loci_lengths){
	
LOOP:	my $random_chr = int(rand($total_chr )); #Returns a random fractional number greater
		# than or equal to 0 and less than the value of EXPR. (EXPR should be positive.)
		my $region_len = $loci_lengths[$i];
		my ($chr, $start, $end);
	
		if ($random_chr < $len_chr1){
			$chr = "chr1";
			$start = $random_chr + 1;
			$end = $start + $region_len - 1;
			if ($end > $len_chr1) {	#$end = $len_chr1
				print STDERR "\n";
				print STDERR join("\t",($i,$random_chr,$chr, $start, $end)),"\n";
				goto LOOP;
			}
		}elsif ( $random_chr < $point_2){
			$chr = "chr2";
			$start = $random_chr - $len_chr1  + 1;
			$end = $start + $region_len - 1;
			if ($end > $len_chr2) {#$end = $len_chr2
				print STDERR "\n";
				print STDERR join("\t",($i,$random_chr,$chr, $start, $end)),"\n";
				goto LOOP;
			}
		}elsif ( $random_chr < $point_3){
			$chr = "chr3";
			$start = $random_chr - $point_2 +1 ;
			$end = $start + $region_len - 1;
			if ($end > $len_chr3) {#$end = $len_chr3
				print STDERR "\n";
				print STDERR join("\t",($i,$random_chr,$chr, $start, $end)),"\n";
				goto LOOP;
			}
		}elsif ( $random_chr < $point_4){
			$chr = "chr4";
			$start = $random_chr - $point_3 +1;
			$end = $start + $region_len - 1;
			if ($end > $len_chr4) {#$end = $len_chr4
				print STDERR "\n";
				print STDERR join("\t",($i,$random_chr,$chr, $start, $end)),"\n";
				goto LOOP;
			}
		}elsif ( $random_chr < $total_chr){
			$chr = "chr5";
			$start = $random_chr - $point_4 +1;
			$end = $start + $region_len - 1;
			if ($end > $len_chr5) {#$end = $len_chr5
				print STDERR "\n";
				print STDERR join("\t",($i,$random_chr,$chr, $start, $end)),"\n";
				goto LOOP;
			}
		}
		else{
			#die "$random_chr >= $total_chr:$!";	
			print STDERR "\n\n\nwrong: $random_chr >= $total_chr \n\n\n\n";
		}
		$ref->[$i] = [$chr, $start, $end];
	}
}


sub read_loci_list{
	my ($list_file, $ref, $ref_length) = @_;
	die "list_file" unless (-e $list_file);
#	print STDERR "reading $list_file...\n";
	open(IN, $list_file) or die "cannot open $list_file";
	my $head = <IN>;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my $len_temp = $a[2] - $a[1] + 1;
		$ref->[$i] = $len_temp;
		$$ref_length += $len_temp;
	}	
#	print STDERR "last_index : $i\n\n";
	close(IN);
}
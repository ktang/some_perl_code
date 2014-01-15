#!/usr/bin/perl -w
# both strand and score
# use promoter regions only, no downstream

use strict;

my $debug = 0;

#my $window_size = 200;
#my $step_size = 100;
#my $min_fold_change = 2;
#my $min_cutoff = 5; #mC_ros3 at least >= $min_cutoff

#my $wig_Col0 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Col0_s_1_methylation.wig";
#my $wig_ros3 = "/Users/tang58/deep_seq_analysis/Ros3_bind_RNA/combined_low_hits_read/IP_uniq_annotated/Ros3_2_s_8_methylation.wig";
#my $wig_Col0 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_col0_C.sgr";
#my $wig_ros3 = "/Users/tang58/Nimblegen_project/Ecker_data/tair8_to_tair9/tair9_mc_rdd_C.sgr";

#if($debug){
#	$wig_Col0 = "/Users/tang58/try/Ros3_binding/wig_files/Col0.wig";
#	$wig_ros3 = "/Users/tang58/try/Ros3_binding/wig_files/Ros3.wig";	
#}

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

#my $usage = "$0 <indir> <outdir> <rdd|ros1|ros3> <fold> <min_number>";
my $usage = "$0 <indir> <outdir> <rdd|ros1|ros3> <min_coverage_for_mut>";
die $usage unless (@ARGV == 4);

my ($indir, $outdir, $label) = @ARGV[0..2];

my $min_coverage = $ARGV[3];

#my ($min_fold_change, $min_cutoff) = @ARGV[3..4];

die "wrong directory\n\n" unless (-e $indir and -e $outdir);

die "rdd|ros1|ros3\n\n" unless ($label eq "rdd" or $label eq "ros1" or $label eq "ros3"
	or $label eq "hsp" or $label eq "phd" or $label eq "zdp" or $label eq "ros2" or $label eq "ros4" );

my $wig_col0 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ1_Col0.wig";

my ($mut1,$mut2);

if($label eq "rdd"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_17_rdd.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_18_rdd.wig";
}elsif($label eq "ros1"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ3_ros1_4.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ4_ros1_4.wig";
}elsif($label eq "ros3"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ7_ros3_2.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ8_ros3_2.wig";
}elsif($label eq "hsp"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ11_hsp.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ12_hsp.wig";
}elsif($label eq "phd"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ12_phd.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ13_phd.wig";
}elsif($label eq "zdp"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_15_zdp.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ_16_zdp.wig";
}elsif($label eq "ros2"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ6_ros2.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ6_ros2.wig";
}elsif($label eq "ros4"){
	$mut1 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ9_ros4.wig";
	$mut2 = "/Users/tang58/deep_seq_analysis/Bisulfite/work_with_Liu/wig_file/JKZ10_ros4.wig";
}
else{
	die "wrong label\n\n";
}

print STDERR join("\n", ($mut1,$mut2, $wig_col0)),"\n\n";

opendir(DIR, $indir) or die "cannot open $indir";

my @files = grep /bed\.txt$/, readdir DIR;

print STDERR join("\n", @files) , "\n\n";

open(COL, $wig_col0) or die "cannot open $wig_col0";
open(M1, $mut1) or die "cannot open $mut1";
open (M2, $mut2) or die "cannot open $mut2";

#my (%cols_p, %cols_m, %muts_p, %muts_m); # distinguish plus strand and minus strand

my (%cols, %muts);

#my ($chr_col, $chr_ros);
my $chr = "";

while(<COL>){
	next if (/^track/);
	if (/variableStep\s+chrom=(\w+)/){
		$chr = lc $1;
		next;
	}
	chomp;
	my ($pos, $val) = split /\t/;
	$cols{$chr}->[$pos] = $val
#	if ($val > 0) { $cols_p{$chr}->[$pos] = $val;}
#	elsif ($val < 0 ){ $cols_m{$chr}->[$pos] = $val;}
#	else {print STDERR "wrong val:$pos\t$val\n"}
#	$cols{$chr}->[$pos] = $val;
}
close(COL);

while(<M1>){
	next if (/^track/);
	if (/variableStep\s+chrom=(\w+)/){
		$chr = lc $1;
		next;
	}
	chomp;
	my ($pos, $val) = split /\t/;

	$muts{$chr}->[$pos] = $val
	
#	if ($val > 0) { $muts_p{$chr}->[$pos] = $val;}
#	elsif ($val < 0 ){ $muts_m{$chr}->[$pos] = $val;}
#	else {print STDERR "wrong val:$pos\t$val\n"}
	
	#$muts{$chr}->[$pos] = $val;
}
close(M1);


while(<M2>){
	next if (/^track/);
	if (/variableStep\s+chrom=(\w+)/){
		$chr = lc $1;
		next;
	}
	chomp;
	my ($pos, $val) = split /\t/;
	$muts{$chr}->[$pos] = $val;
#	if ($val > 0) { $muts_p{$chr}->[$pos] = $val;}
#	elsif ($val < 0 ){ $muts_m{$chr}->[$pos] = $val;}
#	else {print STDERR "wrong val:$pos\t$val\n"}
	#$muts{$chr}->[$pos] = $val;
}
close(M2);

foreach my $file (@files){
	if($file =~ /(\S+)_bed/) {
		my $output = $1."_our_".$label."_minCov".$min_coverage ."_methylation_info.txt";
		if ( (-e "$outdir/$output") or !(-e "$indir/$file")){
			die "input not exists or output exists\n";
		}
		
		open (IN, "$indir/$file") or die "cannot open $file to input:$!";
		open (OUT, ">$outdir/$output") or die "cannot open output $output:$!";
		my @ins = <IN>;
		close(IN);
		
		my $head = $ins[0];
		chomp $head;
		#my $h = "mC_Our".$label."_";
		print OUT join("\t", ($head, "mut_0_1kb", "WT_0_1kb", "Hyper_0_1kb","mut_1_2kb", "WT_1_2kb", "Hyper_1_2kb","mut_0_2kb", "WT_0_2kb", "Hyper_0_2kb", )), "\n";
		
		for(my $i = 1; $i <= $#ins; $i++){
			my ($mC_ros3_1kb, $mC_col0_1kb, $hyper_1kb) = (0, 0, "NULL"); # means 0-1kb
			my ($mC_ros3_2kb, $mC_col0_2kb, $hyper_2kb) = (0, 0, "NULL"); # means 1kb-2kb
			my ($mC_ros3_total, $mC_col0_total, $hyper_total) = (0, 0, "NULL"); # means 0-2kb

			my $this = $ins[$i];
			chomp $this;
			my @a = split "\t", $this;
			my ($chr, $start, $end) = @a[0..2];
			$chr = lc $chr;
			
			if ($a[4] eq "+"){
				
					for (my $j = $start; $j <= $end; $j++){
						if(defined $muts{$chr}->[$j]){$mC_ros3_1kb += abs($muts{$chr}->[$j])}
						if(defined $cols{$chr}->[$j]) {$mC_col0_1kb += abs($cols{$chr}->[$j])}
					}
					
					my ($new_start, $new_end) = ( $start - 1000 , $start - 1);
					for (my $j = $new_start; $j <= $new_end; $j++){
						if(defined $muts{$chr}->[$j]){$mC_ros3_2kb += abs($muts{$chr}->[$j]) }
						if(defined $cols{$chr}->[$j]) {$mC_col0_2kb += abs($cols{$chr}->[$j])}
					}
					
					$mC_ros3_total = $mC_ros3_1kb + $mC_ros3_2kb;
					$mC_col0_total = $mC_col0_1kb + $mC_col0_2kb;
					
					$hyper_1kb = isHyper ($mC_ros3_1kb, $mC_col0_1kb, $min_coverage);
					$hyper_2kb = isHyper ($mC_ros3_2kb, $mC_col0_2kb, $min_coverage);
					$hyper_total = isHyper ($mC_ros3_total, $mC_col0_total , $min_coverage);
					
				
		
				print OUT join("\t", ($this, $mC_ros3_1kb, $mC_col0_1kb, $hyper_1kb, $mC_ros3_2kb, $mC_col0_2kb, $hyper_2kb, $mC_ros3_total, $mC_col0_total, $hyper_total, )), "\n";
			}
			
			elsif ($a[4] eq "-"){
				
				for (my $j = $start; $j <= $end; $j++){
					if(defined $muts{$chr}->[$j]){$mC_ros3_1kb += abs($muts{$chr}->[$j])}
					if(defined $cols{$chr}->[$j]) {$mC_col0_1kb += abs($cols{$chr}->[$j])}
				}
			
				my ($new_start, $new_end) = ( $end + 1, $end + 1000);
				for (my $j = $new_start; $j <= $new_end; $j++){
					if(defined $muts{$chr}->[$j]){$mC_ros3_2kb += abs($muts{$chr}->[$j]) }
					if(defined $cols{$chr}->[$j]) {$mC_col0_2kb += abs($cols{$chr}->[$j])}
				}

				$mC_ros3_total = $mC_ros3_1kb + $mC_ros3_2kb;
				$mC_col0_total = $mC_col0_1kb + $mC_col0_2kb;
					
				$hyper_1kb = isHyper ($mC_ros3_1kb, $mC_col0_1kb, $min_coverage);
				$hyper_2kb = isHyper ($mC_ros3_2kb, $mC_col0_2kb, $min_coverage);
				$hyper_total = isHyper ($mC_ros3_total, $mC_col0_total , $min_coverage);
					

		
				print OUT join("\t", ($this, $mC_ros3_1kb, $mC_col0_1kb, $hyper_1kb, $mC_ros3_2kb, $mC_col0_2kb, $hyper_2kb, $mC_ros3_total, $mC_col0_total, $hyper_total, )), "\n";
			}

			else{
				print STDERR "wrong strand $this\n";
			}
		}		
	}else{
		print STDERR "wrong file $file\n";		
	}
}

exit;

sub isHyper{
	my ($mut, $wt , $cutoff) = @_;

	if($wt == 0){
		if ($mut >= $cutoff){
			return "YES";
		}else{
			return "NO";
		}
	}else{
		if($mut >= $cutoff and ($mut - $wt)/$wt >= 0.2){
			return "YES";
		}
		else{
			return  "NO";
		}
	}			
				
}
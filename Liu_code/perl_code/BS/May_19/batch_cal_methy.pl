#!/usr/bin/perl -w
#
use strict;

my $log = "cal_meth.log";
my %rates = ("Lane1_JKZ1_Col0_meth_db_input.txt"=>0.0014,
             "Lane4_JKZ3_ros1_4_meth_db_input.txt"=>0.0055,
			 "Lane5_JKZ4_ros1_4_meth_db_input.txt"=>0.0055,
			 "Lane6_JKZ12_hsp_meth_db_input.txt"=>0.0059,
			 "Lane7_JKZ6_ros2_meth_db_input.txt"=>0.0049,
			 "Lane8_JKZ7_ros3_2_meth_db_input.txt"=>0.020,
			 "Lane1_JKZ8_ros3_2_meth_db_input.txt"=>0.021,
			 "Lane3_JKZ9_ros4_meth_db_input.txt"=>0.017,
			 "Lane4_JKZ10_ros4_meth_db_input.txt"=>0.020,
			 "Lane5_JKZ11_hsp_meth_db_input.txt"=>0.014,
			 "Lane6_JKZ12_hsp_meth_db_input.txt"=>0.013,
			 "Lane7_JKZ13_phd_meth_db_input.txt"=>0.013,
			 "Lane8_JKZ12_phd_meth_db_input.txt"=>0.023,
			 "JKZ_15_zdp_lane1_meth_db_input.txt"=>0.024,
			 "JKZ_16_zdp_lane2_meth_db_input.txt"=>0.014,
			 "JKZ_17_rdd_lane3_meth_db_input.txt"=>0.012,
			 "JKZ_18_rdd_lane5_meth_db_input.txt"=>0.011);
foreach my $file(sort keys %rates){
	print STDERR "working on file $file ...\n";
	my ($pre, $ext) = split /\./, $file;
	my $outfile = $pre . "_onetype" . "." . $ext;
	my $r = $rates{$file};
	`perl calculate_methy_status_notype.pl $file $r > $outfile 2>>$log`;
}


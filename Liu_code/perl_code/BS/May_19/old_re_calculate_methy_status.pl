#!/usr/bin/perl -w
# recalcualte methylation status using Lister et al. 2009 method
use strict;
use Math::CDF qw(:all);
my $FDR = 0.01;
my $print_values = 0;
my $usage = "$0 <position file> <meth file> <error rate>";
die $usage unless(@ARGV >= 3);
my ($pos_file, $meth_file, $rate) = @ARGV[0..2];

open(PF, $pos_file) or die "Can't open $pos_file: $!";
my @type;
my $n = 0;
my $base_pos = 1;
while(<PF>){
	if(/pos_id/){
		next;
	}
	$n++;
	chomp;
	my @temp = split /\t/;
	if($n == 1){
		$base_pos = $temp[0];
	}
	if($temp[2] =~ /chr\d/){
	    my $index = $temp[0] - $base_pos;
	    $type[$index] = $temp[5];
	}
}
close PF;

my %meth_level; # CG, CHG, or CHH -> [0..n]->[0..n]
my %meth_max;
my %accu_meth;
my %total_meth;
open(MF, $meth_file) or die "Can't open $meth_file: $!";
while(<MF>){
	next if (/SampleID/);
	chomp;
	my ($sample_id, $pos_id, $depth, $num_C, $percent, $isMeth) = split /\t/;
	my $index = $pos_id - $base_pos;
	if(defined $type[$index]){# only consider 5 nuclear chromosomes
		my $t = $type[$index];
		$meth_level{$t}->[$depth]->[$num_C]++;
		$total_meth{$t}->[$depth]++;
		#foreach my $i(0..$num_C){
		#	$accu_meth{$t}->[$depth]->[$i]++;
		#}
		if(!defined $meth_max{$t}){
			$meth_max{$t} = $depth;
		}elsif($meth_max{$t} < $depth){
			$meth_max{$t} = $depth;
		}
	}else{
		#die "type not defined for index $index, pos=$pos_id\n";
	}
}
#undef @type;
close MF;

# calculate accumulated number of positions with <= k unconverted Cs
foreach my $t(keys %meth_max){
	foreach my $i(0..$meth_max{$t}){
		if(!defined $meth_level{$t}->[$i]->[0]){
			$meth_level{$t}->[$i]->[0] = 0;
		}
		$accu_meth{$t}->[$i]->[0] = $meth_level{$t}->[$i]->[0];
		next if($i == 0);
		foreach my $j(1..$i){
			if(!defined $meth_level{$t}->[$i]->[$j]){
				$meth_level{$t}->[$i]->[$j] = 0;
			}
			$accu_meth{$t}->[$i]->[$j] = $accu_meth{$t}->[$i]->[$j-1] +
			     $meth_level{$t}->[$i]->[$j];
		}
	}
}
			
print STDERR "Max depth\n";
foreach my $t(sort keys %meth_max){
	print STDERR $t, " ", $meth_max{$t}, " ";
}
print STDERR "\n";
#die "stop here";


#my $max_depth = 18;
# calculate at each depth, how many positions have k unconverted C
if($print_values){
    print "Number of positions with k unconverted Cs at each depth\n";
}
foreach my $t(sort keys %meth_max){
	my $max_depth = $meth_max{$t};
	if($print_values){
	print $t;
    foreach my $i(0..$max_depth){
		print "\t", $i;
	}
	print "\n";
	}
	foreach my $i(0..$max_depth){
		if($print_values){
		    print $i;
		}
		foreach my $j(0..$i){
		#	if(!defined $meth_level{$t}->[$i]->[$j]){
		#		$meth_level{$t}->[$i]->[$j] = 0;
		#	}
			if($print_values){
			   print "\t", $meth_level{$t}->[$i]->[$j], "|",
			   $accu_meth{$t}->[$i]->[$j];
			}
		}
		if($print_values){
		    print "\n";
		}
	}
}

undef %meth_level;

# find cutoff at each depth level
my %cutoff;
my %num_pos;
if($print_values){
    print "Lister equation values\n";
}
foreach my $t(sort keys %meth_max){
	my $max_depth = $meth_max{$t};
	if($print_values){
	    print $t;
    foreach my $i(0..$max_depth){
		print "\t", $i;
	}
	print "\n";
	}

	foreach my $i(0..$max_depth){
		next if(!defined $total_meth{$t}->[$i]);
		if($print_values){
		    print $i;
		}
		if($i < 2){
			if($print_values){
			    print "\t0";
			}
			$cutoff{$t}->[$i] = 0;
			#next;
		}else{
			$cutoff{$t}->[$i] = 4;
		foreach my $j(1..$i){
			
			my $prob = pbinom($j, $i, $rate);
			if($j > 0){
				$prob = $prob - pbinom($j-1, $i, $rate);
			}
			my ($num_unC, $num_mC) = (0,0);
			$num_unC = $accu_meth{$t}->[$i]->[$j-1];
			if(!defined $num_unC){
				$num_unC = 0;
			}
			$num_mC = $total_meth{$t}->[$i] - $num_unC;
			my $left = $prob * $num_unC;
			my $right = $FDR * $num_mC;
			if($left < $right){
				$cutoff{$t}->[$i] = $j;
				$num_pos{$t}->[$i] = $num_mC;
				last;
			}
			if($print_values){
			    print "\t", $left, "|", $right;
			}
		}
		if($print_values){
		    print "\n";
		}
	}
	}
}

undef %total_meth;
if($print_values){
    print "Cutoff values and number of mC at each depth\n";
}
#my $max_depth = 100;
print STDERR "Total mC position in ", $meth_file, "\n";
foreach my $t(sort keys %meth_max){
	my $max_depth = $meth_max{$t};
	#if($print_values){
	    #print STDERR $t, "\n";
	#}
    foreach my $i(0..$max_depth){
		if($print_values){
		print $i;
		if($i != $max_depth){
			print "\t";
		}
		}
	}
	if($print_values){
	    print "\n";
	}
	if($print_values){
	foreach my $i(0..$max_depth){
		
		if(defined $cutoff{$t}->[$i]){
			print $cutoff{$t}->[$i];
		}else{
			print "0";
		}
		if($i != $max_depth){
			print "\t";
		}
	}
	print "\n";
	}
	my $total = 0;
    foreach my $i(0..$max_depth){
		if(defined $num_pos{$t}->[$i]){
			if($print_values){
			    print $num_pos{$t}->[$i];
			}
			$total += $num_pos{$t}->[$i];
		}else{
			if($print_values){
			    print "0";
			}
		}
		if($i != $max_depth){
			if($print_values){
			    print "\t";
			}
		}
	}
	if($print_values){
	    print "\n";
	}
	print STDERR "$t: $total\n";

}
undef %meth_max;
undef %num_pos;
open(MF, $meth_file) or die "Can't open $meth_file: $!";
#seek(MF, 0, 0) or die "Can't seek to beginning of file: $!";
while(<MF>){
	if (/SampleID/){
		print;
		next;
	}
	chomp;
	my ($sample_id, $pos_id, $depth, $num_C, $percent, $isMeth) = split /\t/;
	my $t = $type[$pos_id - $base_pos];
	if(defined $t){
		$isMeth = 0;
		if($depth >= 2 && $num_C >= 1 && $num_C >= $cutoff{$t}->[$depth]){
			$isMeth = 1;
		}
		print join("\t", ($sample_id, $pos_id, $depth, $num_C, $percent, 
		    $isMeth)), "\n";
	}
}
close MF;

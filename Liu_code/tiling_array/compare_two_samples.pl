#!/usr/bin/perl -w
# compare probe intensities in two samples (a mutant vs control) and identify regions that show higher or
# lower expression in mutant
use strict;
use Statistics::Descriptive;

my $min_sig_probes = 5;
my $max_dis_btw_probes = 100; #bp
my $probe_length = 25;
my $max_hypo_diff = -1; #log2 difference
my $min_hyper_diff = 1; #log2 difference

my $usage = "$0 <probe intensity file> <col index 1> <col index 2> <1 for hypo, 0 for hyper>";
die $usage unless (@ARGV >= 4);
my ($valFile, $index1, $index2, $isHypo) = @ARGV[0..3];

my (%orderedPos, %vals);
open(VF, $valFile) or die "Can't open $valFile:$!";
my @header;
while(<VF>){
	chomp;
	if(/Start/){
		@header = split /\t/;
		next;
	}
	my @temp = split /\t/;
	my ($chr, $start) = @temp[0..1]; 
	push @{$orderedPos{$chr}}, $start;
	$vals{$chr}->{$start} = [$temp[$index1], $temp[$index2]];	
}

print join("\t", ("Chr", "Start", "End", "num_sig_probes", "num_total_probes", $header[$index1], $header[$index2], "Log2_Val_diff")), "\n";

my (@accum1, @accum2, @inter1, @inter2);
my ($num_sig, $num_total, $last_start, $last_end, $last_chr) = (0, 0, 0, 0, "NONE");
foreach my $chr(sort keys %orderedPos){
	if(($chr ne $last_chr) && ($last_chr ne "NONE")){ #print last region on the last chromosome if necessary
		if($num_sig >= $min_sig_probes){
			if($num_total != scalar(@accum1)){
				die "Number of total probes does not equal to the array size";
			}
			print_vals($last_chr, $last_start, $last_end, $num_sig, $num_total, \@accum1, \@accum2);
		}
		@accum1 = (); @inter1 = ();
		@accum2 = (); @inter2 = ();
		($num_sig, $num_total, $last_start, $last_end, $last_chr) = (0, 0, 0, 0, $chr); 
	}
			
	my @pos = @{$orderedPos{$chr}};
	foreach my $i(0..$#pos){
		my $diff = $vals{$chr}->{$pos[$i]}->[1] - $vals{$chr}->{$pos[$i]}->[0];
		if($isHypo){# hypo
			if($diff <= $max_hypo_diff){# a sig probe
				if($last_end == 0){#first sig probe
					$last_start = $pos[$i];
				}else{
					my $dist = $pos[$i] - $last_end;
					if($dist <= $max_dis_btw_probes){
						if(scalar(@inter1) >= 1){
							push @accum1, @inter1;
							push @accum2, @inter2;
							$num_total += scalar(@inter1);
							@inter1 = ();
							@inter2 = ();
						}
					}else{ #write out and start over
					    if($num_sig >= $min_sig_probes){
			                print_vals($chr, $last_start, $last_end, $num_sig, $num_total, \@accum1, \@accum2);
		                }
		                @accum1 = (); @inter1 = ();
		                @accum2 = (); @inter2 = ();
						$num_sig = 0; $num_total = 0;
						$last_start = $pos[$i];
					}
				}
				$last_end = $pos[$i]+$probe_length-1;	
				push @accum1, $vals{$chr}->{$pos[$i]}->[0];
				push @accum2, $vals{$chr}->{$pos[$i]}->[1];
				$num_sig++; $num_total++;
			}else{ #not a sig probe
			    if($num_sig > 0){
					push @inter1, $vals{$chr}->{$pos[$i]}->[0];
					push @inter2, $vals{$chr}->{$pos[$i]}->[1];
				}
			}

		}else{# hyper
		    	if($diff >= $min_hyper_diff){# a sig probe
				if($last_end == 0){#first sig probe
					$last_start = $pos[$i];
				}else{
					my $dist = $pos[$i] - $last_end;
					if($dist <= $max_dis_btw_probes){
						if(scalar(@inter1) >= 1){
							push @accum1, @inter1;
							push @accum2, @inter2;
							$num_total += scalar(@inter1);
							@inter1 = ();
							@inter2 = ();
						}
					}else{ #write out and start over
					    if($num_sig >= $min_sig_probes){
			                print_vals($chr, $last_start, $last_end, $num_sig, $num_total, \@accum1, \@accum2);
		                }
		                @accum1 = (); @inter1 = ();
		                @accum2 = (); @inter2 = ();
						$num_sig = 0; $num_total = 0;
						$last_start = $pos[$i];
					}
				}
				$last_end = $pos[$i]+$probe_length-1;	
				push @accum1, $vals{$chr}->{$pos[$i]}->[0];
				push @accum2, $vals{$chr}->{$pos[$i]}->[1];
				$num_sig++;$num_total++;
			}else{ #not a sig probe
			    if($num_sig > 0){
					push @inter1, $vals{$chr}->{$pos[$i]}->[0];
					push @inter2, $vals{$chr}->{$pos[$i]}->[1];
				}
			}

		}
}			
}
if($num_sig >= $min_sig_probes){
			if($num_total != scalar(@accum1)){
				die "Number of total probes does not equal to the array size";
			}
			print_vals($last_chr, $last_start, $last_end, $num_sig, $num_total, \@accum1, \@accum2);
		}

#sub print_vals{
#	my ($chr, $start, $end, $num_s, $num_t, $ref1, $ref2) = @_;
#	my @vals1 = @$ref1;
#	my @vals2 = @$ref2;
#   my $stat1 = Statistics::Descriptive::Full->new();
#	$stat1->add_data(@vals1);
#	my $trimmed_mean1 = $stat1->trimmed_mean(0.2);
#	my $stat2 = Statistics::Descriptive::Full->new();
#	$stat2->add_data(@vals2);
#	my $trimmed_mean2 = $stat2->trimmed_mean(0.2);
#	my $diff = $trimmed_mean2 - $trimmed_mean1;
#	print join("\t", ($chr, $start, $end, $num_s, $num_t, $trimmed_mean1, $trimmed_mean2, $diff)), "\n";
#}
sub print_vals{
	my ($chr, $start, $end, $num_s, $num_t, $ref1, $ref2) = @_;
	my @vals1 = @$ref1;
	my @vals2 = @$ref2;
	my $stat1 = Statistics::Descriptive::Full->new();
	$stat1->add_data(@vals1);
	my $trimmed_mean1 = $stat1->trimmed_mean(0.2);
	my $stat2 = Statistics::Descriptive::Full->new();
	$stat2->add_data(@vals2);
    my $trimmed_mean2 = $stat2->trimmed_mean(0.2);

    my $diff_stat = Statistics::Descriptive::Full->new();
	foreach my $i(0..$#vals1){
	    $diff_stat->add_data($vals2[$i] - $vals1[$i]);
	}
	
	my $diff = $trimmed_mean2 - $trimmed_mean1;
	my $mean_diff = $diff_stat->mean();
	if(($isHypo && ($mean_diff <= $max_hypo_diff)) || (!$isHypo && ($mean_diff >= $min_hyper_diff))){ 
	    print join("\t", ($chr, $start, $end, $num_s, $num_t, $trimmed_mean1, $trimmed_mean2, $diff)), "\n";
	}
}


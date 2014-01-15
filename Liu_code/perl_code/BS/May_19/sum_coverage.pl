#!/usr/bin/perl -w
# summarize coverage info
use strict;

my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);

print join("\t", ("Sample", "PercentCoverage>=2Watson", "PercentCoverage>=2Crick", "AverageDepthWatson", "AverageDepthCrick")), "\n";

foreach my $dir(@ARGV){
    opendir(CDIR, $dir);
    my @dirs = sort readdir CDIR;
    foreach my $d(@dirs){
		my $this_dir = "$dir/$d";
        if((-d $this_dir) && ($this_dir =~ /Lane/i)){
		    print STDERR "calculating coverage under $this_dir\n";
		    opendir(SDIR, $this_dir);
		    my @files = grep {/coverage_stats/} readdir SDIR;
			if(@files >= 2){
				die "more than two coverage files";
			}
		    foreach my $file(@files){
				$file = $this_dir . "/" . $file;
				open(IN, $file) or die "Can't open $file: $!";
	            my ($in_cover, $in_depth) = (0, 0);
				my (@gt_two, %depth);
				while(<IN>){
					if(/Coverage >= 2/){
						$in_cover = 1;
						next;
					}
					if(/Depth coverage/){
						$in_cover = 0;
						$in_depth = 1;
						next;
					}
					if($in_cover){
						if(/Whole genome\s+(\d+\.\d+)/){
							push @gt_two, $1;
						}
					}
					if($in_depth){
						if(/(chr\d)\s+(\d+\.\d+)/){
							if(defined $depth{$1}){
								$depth{$1}->[1] = $2;
							}else{
								$depth{$1}->[0] = $2;
							}
						}
					}
				}
					my ($t_watson, $t_crick, $total_len) = (0,0, 0);
					foreach my $chr(keys %depth){
						$t_watson += $depth{$chr}->[0] * $chr_len{$chr};
						$t_crick += $depth{$chr}->[1] * $chr_len{$chr};
						$total_len += $chr_len{$chr};
					}
					my $avg_watson = $t_watson / $total_len;
					my $avg_crick = $t_crick / $total_len;
					print join("\t", ($d, @gt_two, $avg_watson, $avg_crick)),
					           "\n";
				
			}
	 
	}
}
}

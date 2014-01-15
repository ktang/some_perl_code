#!/usr/bin/perl -w

#Dec 6, 2013
#modified from
#/Users/tang58/Kai_BS/myown/cytosines_coverage/cal_cytosines_coverage_for_isMeth_output_v0.0.pl
# input isMeth file
# output two files for meth level notEcker and Ecker
#column: num_with_dep CG CHG CHH C
#row: chr1-5,C,M,... +
#     ...............-
#chr1-5
# require depth cutoff (usually 4)

########################################################################################
# given acgt-count output, calculate the depth of coverage at each cytosine position
# and percentage of total cytosines that was covered by at least one or two reads.

#chr1    33      33      CHG:13  0.153846        -
#chr1    79      79      CHH:27  0.148148        -
# 0		  1		 2			3		4			 5

#v0.1 output depth file

use strict;
use File::Spec;

#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502, "chrM"=>366924, "chrC"=>154478, "pUC19"=>2686);


my $debug = 0;
my $constant = 2000000;
if ($debug) {$constant = 300 ;}

if($debug){
	print  STDERR "debug = $debug\n\n";
}

#my $usage = "$0 \n <indir> <infile> <outdir> <pre> \n\n";
#die $usage unless (@ARGV == 4);
my $dep_cutoff = 4;

my $usage = "$0 \n <isMeth_file> <outdir> <outpre> [<dep_cutoff>]\n\n";
die $usage unless (@ARGV >= 3);

my $input = shift or die "input";

#my $indir = shift or die "indir";
#my $infile = shift or die "infile";
my $outdir = shift or die "outdir";
my $pre = shift or die "pre";

if (@ARGV>=1) {
	$dep_cutoff = shift or die "cutoff";
	print STDERR "default cutoff is 4, input cutoff is $dep_cutoff\n\n";
}


#die "wrong indir" unless (-d $indir);
die "wrong outdir" unless (-d $outdir);

#my $in_forw = File::Spec->catfile($indir, $pre . "_forw.txt");
#my $in_rev  = File::Spec->catfile($indir, $pre . "_rev.txt");
#my $input = File::Spec->catfile($indir, $infile);
die unless (-e $input);



if($debug){
	print  STDERR "input files:\n";
#	print  STDERR join("\n", ($in_forw, $in_rev)), "\n\n";
	print  join $input, "\n\n";
}
#die "wrong input" unless (-e $in_forw and -e $in_rev );

#my $output = File::Spec->catfile($outdir,$pre. "_cytosines_coverage_info_DepthCutoff.txt");
my $output_E  = File::Spec->catfile($outdir,$pre. "_meth_level_Ecker_Dep" .$dep_cutoff.   ".txt");
my $output_nE = File::Spec->catfile($outdir,$pre. "_meth_level_notEcker_Dep" .$dep_cutoff.".txt");
#die "$output exists" if (-e $output);


if($debug){
	#print STDERR "output:\n$output\n";
}

my @C_type = ( "C", "CG", "CHG", "CHH");

my @mmC_type = map { "mm_" . $_} @C_type; # mmC = sum of per /C number
my @wmC_type = map { "wm_" . $_} @C_type;
my @num_dep_type = map {$_ ."_num_with_dep". $dep_cutoff }  @C_type;
#"num_with_dep". $dep_cutoff

die "$output_E exists" if (-e $output_E);
open(OUT_E, ">>$output_E" ) or die;
print OUT_E join("\t", ("chr", "strand", @num_dep_type ,  @wmC_type, @mmC_type)), "\n";
close OUT_E;
 

die "$output_nE exists" if (-e $output_nE);
open(OUT_NE, ">>$output_nE" ) or die;
print OUT_NE join("\t", ("chr", "strand", @num_dep_type,  @wmC_type, @mmC_type)), "\n";
close OUT_NE;
#my $depth_file = File::Spec->catfile($outdir,$pre. "_depth_distribution_DepthCutoff.txt" );
#die "depth_file exists" if (-e $depth_file);

#my $mC_depth_file = File::Spec->catfile($outdir,$pre. "_depth_mC_DepthCutoff.txt");
#die if (-e $mC_depth_file);

if($debug){
	#print STDERR $depth_file, "\n\n";
}
if($debug){
	print STDERR "\nOK\n\n";
	exit;
}

#die "die debug\n\n" if ($debug);

#open (OUT, ">>$output") or die;
#open (DEP, ">>$depth_file") or die;
#print DEP join("\t", ("depth", "number")), "\n";

my %num_with_dep;
my %deps_record; #deps{$chr}->{strand}->{type} += 

my %mC_Ecker;
my %mC_notEcker;

my %sum_per_Ecker;
my %sum_per_notEcker;


open(IN, $input) or die;
# 0	1	2	  3	4	5	6		7
#chr     pos     strand  type    num_C   depth   percentage      isMeth
#chr1    334     +       CHH     0       20      0       0
#chr1    336     +       CHH     0       20      0       0

my $head = <IN>;

while(<IN>){
	chomp;
	my @a = split "\t";
	
	my $depth = $a[5];
	next unless ( $depth >= $dep_cutoff );
	
	my $chr    = $a[0];
	my $strand = $a[2];
	my $type   = $a[3];
	my $C_num_original = $a[4];
	my $per    = $a[6];
	my $isMeth = $a[7];
	
	$mC_notEcker{$chr}->{$strand}->{$type} += $C_num_original;
	$mC_notEcker{$chr}->{$strand}->{C} += $C_num_original;
	
	$sum_per_notEcker {$chr}->{$strand}->{$type} +=  $per;
	$sum_per_notEcker {$chr}->{$strand}->{C} +=  $per;
	
	if ( $isMeth == 1) {
		$mC_Ecker{$chr}->{$strand}->{$type} += $C_num_original;
		$mC_Ecker{$chr}->{$strand}->{C} += $C_num_original;
		
		$sum_per_Ecker{$chr}->{$strand}->{$type} +=  $per;
		$sum_per_Ecker{$chr}->{$strand}->{C} +=  $per;
	}
	
	
	$num_with_dep{$chr}->{$strand}->{$type}++;
	$num_with_dep{$chr}->{$strand}->{C}++;
	
	$deps_record{$chr}->{$strand}->{$type} += $depth;
	$deps_record{$chr}->{$strand}->{C} += $depth;

	
}

close (IN);

output_file($output_E,  \%num_with_dep, \%deps_record, \%mC_Ecker,    \%sum_per_Ecker);
output_file($output_nE, \%num_with_dep, \%deps_record, \%mC_notEcker, \%sum_per_notEcker);


exit;

#output_file($output_E,  \%num_with_dep, \%deps_record, \%mC_Ecker,    \%sum_per_Ecker);

sub output_file{
	my ( $file, $num_ref, $dep_ref, $mC_ref, $sum_per_ref ) = @_;
	my %chrs = ("chr1" => 1,
		    "chr2" => 1,
		    "chr3" => 1,
		    "chr4" => 1,
		    "chr5" => 1,
		    );
	my @types = ("C", "CG", "CHG", "CHH");

	open(OUT, ">>$file") or die;
	
	my ( %num_all, %dep_all, %mC_all, %sum_per_all )  = (0) x 4;
	foreach my $chr (sort keys %{$num_ref}){
		foreach my $strand (sort keys %{ $num_ref->{$chr} }){
			# print first two cols
			print OUT join("\t", ( $chr, $strand) ), "\t";
			
			my @nums = ();
			my @levels = ();
			my @mm_per = ();
			
			for my $type (@types){
				my ( $num, $dep, $mC, $sum_per )  = (0) x 4;
				if(defined $num_ref->{$chr}->{$strand}->{$type}){
					$num		= $num_ref->{$chr}->{$strand}->{$type}
				}
				 if(defined $dep_ref->{$chr}->{$strand}->{$type}){
					$dep		= $dep_ref->{$chr}->{$strand}->{$type}
				}
				if(defined $mC_ref->{$chr}->{$strand}->{$type}){
					$mC		= $mC_ref->{$chr}->{$strand}->{$type}
				}
				if(defined $sum_per_ref->{$chr}->{$strand}->{$type}){
					$sum_per	= $sum_per_ref->{$chr}->{$strand}->{$type}
				}
				push @nums, $num;
				if ($dep != 0) {
					my $level    = sprintf ("%.3f", 100 * $mC / $dep) . "%";
					my $mm_level = sprintf ("%.3f", 100 * $sum_per / $num) . "%";
					
					push @levels, "$mC/$dep=$level";
					push @mm_per, "$sum_per/$num=$mm_level";
					
				}else{
					push @levels, "NA";
					push @mm_per, "NA"
				}
				if (defined $chrs{$chr}) {
					$num_all{$type}		+= $num;
					$dep_all{$type}		+= $dep;
					$mC_all{$type}		+= $mC;
					$sum_per_all{$type}	+= $sum_per;
				}
				
			}
			print OUT join("\t", (@nums, @levels, @mm_per)),"\n";

		}
	}
	
	#cal for 5 chrs
	print OUT join("\t", ("chr1-5", "+/-" )), "\t";
	my @nums_all = ();
	my @levels_all = ();
	my @mm_per_all = ();
	for my $type (@types){
		if (defined $dep_all{$type}) {
			
			push @nums_all, $num_all{$type}	;
			
			my $level    = sprintf ("%.3f", 100 * $mC_all{$type}/$dep_all{$type}) . "%";
			my $mm_level = sprintf ("%.3f", 100 * $sum_per_all{$type}/$num_all {$type}) . "%";
			
			push @levels_all, "$mC_all{$type}/$dep_all{$type}=$level";
			push @mm_per_all, "$sum_per_all{$type}/$num_all{$type}=$mm_level";
		}else{
			push @nums_all,  "NA";
			push @levels_all, "NA";
			push @mm_per_all, "NA"
		}
		
	}
	print OUT join("\t", (@nums_all, @levels_all, @mm_per_all)),"\n";
	
	close OUT;
}

sub output_mC_depth{
	my ($file, $CG_ref, $CHG_ref, $CHH_ref, $max) = @_;
	
	die if(-e $file);
	open(DEP, ">>$file") or die;
	print DEP join("\t", ("depth", "mC", "CG", "CHG", "CHH")), "\n";
	for my $i(0..$max){
		for my $j(0..$i){
			my ($num_cg, $num_chg, $num_chh) = (0, 0, 0);
			if(defined $CG_ref ->[$i]->[$j])  {$num_cg = $CG_ref->[$i]->[$j]}
			if(defined $CHG_ref->[$i]->[$j])  {$num_chg = $CHG_ref->[$i]->[$j]}
			if(defined $CHH_ref->[$i]->[$j])  {$num_chh = $CHH_ref->[$i]->[$j]}
			if( $num_cg + $num_chg + $num_chh  > 0){
				print DEP join("\t", ( $i, $j, $num_cg, $num_chg, $num_chh )), "\n";
			}
		}
	}
	
	close(DEP);
}



sub round {
    my ($number) = shift;
    #return int($number + .5);
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}
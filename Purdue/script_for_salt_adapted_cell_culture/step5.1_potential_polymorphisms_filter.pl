#!/usr/bin/perl -w

#step6_categorize_variation_filter_and_count_for_each_lib.pl
# Purpose: analyze file Jan8_var_db_pos_14240loci_in_19bam_pileup_phred20_dep8_noN_categorize.txt
# which have format like
# 0	 1	2	3
#chr     pos     ref     N01_S27_XZ9_Col0_E      N02_S26_XZ2_Col0_0mM	N03_B1_real_Col0_0mM
#1       4290    C       DEP=12;TYPE=DEL;ALT=DEL;PER=8.33        DEP=10;TYPE=INS;ALT=INS;PER=10.00       DEP=6;TYPE=NA

#DEP=0   DEP=7;TYPE=NA   DEP=19;TYPE=DEL;ALT=DEL;PER=5.26

# output format shold be like
#


# Transitions: A<=>G C<=>T
# A=>G
# G=>A
# C=>T
# T=>C

#AT=>GC
#GC=>AT


# Transversions: 
# A=>C A=>T C=>G T=>G
# C=>A T=>A G=>C G=>T

#AT=>CG
#AT=>TA
#GC=>TA
#GC=>CG

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;

my @Transitions   = ( "A=>G", "T=>C", "G=>A", "C=>T" ) ;
my @Transversions = ( "A=>C", "T=>G", "A=>T", "T=>A", "G=>T", "C=>A" , "G=>C", "C=>G") ;

my @bidir_a = (  "AT=>GC", "GC=>AT",   "AT=>CG", "AT=>TA", "GC=>TA", "GC=>CG"  ) ;

my %bidir_h = ( "A=>G" => "AT=>GC",
	        "T=>C" => "AT=>GC",
		"G=>A" => "GC=>AT",
		"C=>T" => "GC=>AT",
		
		"A=>C" => "AT=>CG",
		"T=>G" => "AT=>CG",
		"A=>T" => "AT=>TA",
		"T=>A" => "AT=>TA",
		"G=>T" => "GC=>TA",
		"C=>A" => "GC=>TA",
		"G=>C" => "GC=>CG",
		"C=>G" => "GC=>CG"
);

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#my $min_base_qual = 20;

my $min_dep = 8;
my $min_per = 30;

my $usage = "$0 \n <input> <output_list_no_polymor> <filted_lines> <sample_num> \n\n";
die $usage unless(@ARGV == 4);

my $input = shift or die;
my $output = shift or die;
my $filted_file = shift or die;
#my $outpre = shift or die;
my $sample_num  = shift or die;
#my $phred_33_64 = shift or die;
 
#die "Phred_33_64" unless ($phred_33_64 == 33 or $phred_33_64 == 64);
#my $output = $outpre . "_dep" . $min_dep . "_per" . $min_per . ".txt";
die unless (-e $input);
die if( -e $output);
die if( -e $filted_file);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";
open( POL, ">$filted_file") or die;
my $head = <IN>;
print OUT $head;
print POL $head;


while (<IN>) {
	chomp;
	my @a = split "\t";
	
#	print OUT join ("\t", @a[0..2]), "\t";
	#my $ref_base = $a[2];
	my $flag = 0; # if flag = 0; DO NOT output the line; elsif flag = 1 , output it
	my $benchmark = undef;

	for my $n ( 1..$sample_num ){
		my $this = $a[$n + 2];
		my @paires = split ";", $this;
	#DEP=14;TYPE=DEL;ALT=DEL;PER=7.14
	#DEP=17;TYPE=SNP,DEL,INS;ALT=A,DEL,INS;PER=5.88,5.88,11.76
		
	
		if ( @paires==4 ) {
			if( $this =~ /DEP=(\d+);TYPE=(\S+);ALT=(\S+);PER=(\S+)$/ ){
				my ($dep, $type, $alt , $per ) = ($1, $2, $3, $4);
				if ( $dep >= $min_dep) {
					my @types = split ",", $type;
					my @alts = split ",", $alt;
					my @pers = split ",", $per;
					my $ind_of_maximum = Kai_Module::get_ind_of_maximum(\@pers);
					my ($type_max, $alt_max , $per_max ) = ( $types[$ind_of_maximum], $alts[$ind_of_maximum], $pers[$ind_of_maximum]) ;
					
					if ($per_max < $min_per ) {
						$flag = 1;
						last;		
					}else{
						unless ( defined $benchmark ) {
							$benchmark  = $alt_max;#code
						}
						if ($alt_max ne  $benchmark ) {
							$flag = 1;
							last;
						}
						
						
					}
					
				}
			}else{
				die $_;
			}
		}elsif(  @paires==2 ){
			$flag = 1;
			last;
		}elsif(  @paires==1 ){
			next;
		}else{
			die $this;
		}
	}
	
	if ($flag) {
		print OUT $_, "\n";#code
	}else{
		print POL $_, "\n";#code
	}
}
close(IN);
close(OUT);

exit;
#!/usr/bin/perl -w

#v0.0 uncompleted Mar10, 2013

# v0.1
# foreach locus categorize it

#step3.2_categorize_mutation_SNP_to_find_salt_induced_mutation.pl
#only mutated in all three 150 but not 0.
#Also to check whether there is mutated in all 0 but not 150.

#input 
#/Users/tang58/misc/Zhu_Xiaohong/downstream/dep5/step1_all_common_and_other/not_all_common_3404loci.txt
# no head
# chr1    88879   C   ||    17      NONE    NONE    0   ||    15      C=>G    40.00   6    ||   18      C=>G    100.00  18    ||  8       C=>G    100.00  8  ||  16      C=>G    100.00  16    ||  16      C=>G    100.00  16
#  0	   1	  2	  3-6 WT_0	max_num	  	  7-10 WT_150	  	   	       11-14	 ddc_0 	  	       15-18	ddc_150                   19-22 nrpe_0				23-26
#...,.,.,..,,.,,.^K.     ,,...GG.GGG,g,, GGGgGgggGGGGGGGGGg      gGGGGGGG        ggGGGGGgGgGgGGGG        ggggGGggGGgGggGg


#~/misc/Zhu_Xiaohong/downstream/dep5/step1_all_common_and_other_13:23:01_N=522$
#less not_all_common_3404loci.txt | perl ../../src/step3.2_categorize_mutation_SNP_to_find_salt_induced_mutation.pl 2> /dev/null | less
					
use strict;

my $debug = 1;
my %records;

my %refs;

while (<>){
	chomp;
	my @a = split "\t";
	my $indel = 0;
	my ( $chr, $pos, $ref ) = @a[0..2];
	
	my @sample_0   = ( $a[4], $a[12], $a[20] );
	my @sample_150 = ( $a[8], $a[16], $a[24]);
	
	my $flag_0   = decide ( \@sample_0   );
	my $flag_150 = decide ( \@sample_150 );
	# all_mut = 1;
	#all none = -1;
	# else   = 0
	if( $flag_0 == 1 and $flag_150 == -1){
		print STDERR "reverse:", $_, "\n";
	}
	
	if( $flag_0 == -1 and $flag_150 == 1){
		print $_, "\n";
	}
	
}

exit;

# decide ( \@sample_0   );
sub decide{
	my ( $ref )=  @_;
	my @temp = @{$ref};
	
	my $flag_none = 1;
	
	for my $i(0..$#temp){
		if( $temp[$i] ne "NONE" ){
			$flag_none = 0;
		}
	}
	
	my $flag_mut = 1;
	for my $i(0..$#temp){
	#	if( $temp[$i] !~  /=>/ ){
		if( $temp[$i] eq "NONE" ){
			$flag_mut = 0;
		}
	}
	if( $flag_none and  !$flag_mut){
		return -1;
	}elsif( !$flag_none  and $flag_mut ){
		return 1;
	}elsif( $flag_none  and $flag_mut ){
		print STDERR join("\t", @temp), "\n";
		die;
	}else{
		return 0;
	}
}
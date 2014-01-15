#!/usr/bin/perl -w
use strict;
use File::Spec;

#extract_sRNA_db_for_coordinate_list_v0.0.pl
#v0.0
#simply print all list

#v0.1
# for specific length sRNA fold change

my $debug = 0;

# 8/2 = 4
# 8:dividend
#2: divisor
#my $usage = "$0 \n <input_db_mut> <in_db_wt>  <outdir> <outpre> <exp_cutoff>  <fold_change> \n\n";
#die $usage unless(@ARGV == 6);


my $usage = "$0 \n <input_coordinate_data> <db_for_extraction> <label> <output>  \n\n";
die $usage unless(@ARGV == 4);



my $coor_list	= shift or die ;
my $db_file	= shift or die;

my $label 	= shift or die;
my $output 	= shift or die;

die unless ( -e $coor_list);
die unless ( -e $db_file);
die if	   ( -e $output);

my %records;

record( $coor_list, \%records);
die if	   ( -e $output);
open(OUT, ">>$output") or die;

open(DB, $db_file) or die;

while(<DB>){
	print OUT $_ if(/^#/);
	chomp;
	my @a = split "\t";
	if($a[0] eq "coordinate"  ){
		my @b = map { $_ . "_" . $label } @a;
		print OUT join ("\t", @b), "\n";
	}elsif(defined $records{$a[0]}){
		print OUT join( "\t", @a ), "\n";

	}else{
		next;
	}
}
close IN;
close DB;

exit;

#record( $coor_list, \%records);
#  #WT     #reads_for_normalization:5617030        normalization_factor:1.78030026544277
#  #mut    #reads_for_normalization:4608991        normalization_factor:2.1696722775115
#  coordinate      total_exp       fold_change  
sub record{
	my ($file, $ref) = @_;
	die unless (-e $file);
	
	open (IN, $file) or die;
	while(<IN>){
		next if(/^#/);
		chomp;
		my @a = split "\t";
		next if($a[0] eq "coordinate" );
		$ref->{$a[0]} = 1;
	}
	close IN;
}
#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: add_overlap_info_coordinate_format.pl
#
#        USAGE: ./add_overlap_info_coordinate_format.pl  
#
#  DESCRIPTION: input inupt list and a db file, add a col in input to indicate 
#               whether the inverval is exists in db file
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: --- 
#         AUTHOR: Kai Tang (KT), tang58@purdue.edu 
#         ORGANIZATION: Purdue University 
#         VERSION: 1.0
#      CREATED: 07/09/2013 17:21:25
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $debug = 1;

print STDERR "key is X:a-b\n\n";

my $usage = "$0 \n <db_file_used_as_benchmark> <db_label> <input_list> <output _based_on_input_list>\n\n";
die $usage unless (@ARGV == 4);

my $db_file = shift or die;
my $label = shift or die;
my $input   = shift or die;
my $output  = shift or die;

die unless ( -e $db_file);
die unless ( -e $input);
die if ( -e $output);

my %rec;
record_db($db_file,\%rec);

my $n = 0;
open(IN, $input) or die;
open (OUT, ">>$output") or die;
while(<IN>){ 
    if(/^#/){
        print OUT ;
        next;
    }
    chomp;
    if(/^coor/){
        print OUT join("\t", ($_, "overlap_" . $label)), "\n";
        next;
    }
    my @a = split "\t";
    if(defined $rec{$a[0]}){
        $n++;
        print OUT join("\t", ($_, "Y" )), "\n";
    }else{
        print OUT join("\t", ($_, "N" )), "\n";

    }
}
close IN;
close OUT;
print STDERR "overlap_num: $n\n\n";
exit; 

sub record_db{
    my ($file, $ref) = @_;
    die unless (-e $file);
    open (IN, $file) or die;
    while (<IN>){
        chomp;
        my @a = split "\t";
        $ref->{$a[0]} = 1;
    }
    close IN;
}


##WT     #reads_for_normalization:5617030        normalization_factor:1.78030026544277 
###mut    #reads_for_normalization:4608991        normalization_factor:2.1696722775115 
##coordinate      total_exp  

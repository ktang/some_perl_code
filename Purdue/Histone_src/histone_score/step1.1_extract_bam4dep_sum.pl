#!/usr/bin/perl -w

# step1.1_extract_bam4dep_sum.pl
# input bed like file
# add two col , one is sum of dep in intreval and the other is sum/length
# diff from step1 script

use strict;
use File::Spec;

my $debug = 0;
my $print_debug = 0;

if ( !$debug ) {
	$print_debug = 0;
}


#my $half_interval_length = 2000;
#my $half_flanking_bin_num = 100; #40;
#my $bin_size = 50;


if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <input> <output>\n\n";
#die $usage unless(@ARGV == 2);

#my $usage = "$0 \n <input_list_DMR> <bam_dir> <outdir> <output_pre>\n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <input_list_DMR> <bam_dir> <output>\n\n";
die $usage unless(@ARGV == 3);
# bam file: /Volumes/Macintosh_HD_2/Histone_modification/GSE28398/display_simple_bam_ln 

my $input	= shift or die;
my $bam_dir	= shift or die;
my $output	= shift or die;
#my $outdir	= shift or die;
#my $outpre	= shift or die;

#my $output	= shift or die;
#my $output = File::Spec->catfile($outdir, $outpre . "_half" . $half_interval_length . "bp" . "_bin" . $bin_size ."bp.txt");
#my $output = File::Spec->catfile($outdir, $outpre . "_halfBinN" . $half_flanking_bin_num  . "_bin" . $bin_size ."bp.txt");

die unless (-e $input);
die if( -e $output);
die unless (-d $bam_dir);
die if(-e $output);

###########
# get bam file and label
#############
my @bam_files_simple;
my @bam_files_full;

opendir (DIR, $bam_dir);
@bam_files_simple = grep /\.bam$/, readdir DIR;
closedir DIR;

@bam_files_full = map { File::Spec->catfile($bam_dir, $_ )} @bam_files_simple;

my $file_num = $#bam_files_simple + 1;

foreach my $tmp ( @bam_files_full ){
	die $tmp unless (-e $tmp);
}

#print OUT "#";
#print OUT join("\t",@bam_files_simple ), "\n";

my @label = get_bam_label(@bam_files_simple);

if ($debug) {
	print STDERR join("\n", @bam_files_full ),"\n\n";
	print STDERR join("\n", @bam_files_simple ),"\n\n";
	print STDERR join("\n", @label ),"\n\n";
	
	exit;#code
}


###########
# 1, read_list_and_get_interval
##############
open(OUT, ">$output") or die "cannot open $output: $!";

open(IN, $input) or die "input";
my $head = <IN>;
chomp $head;
print OUT $head, "\t";

print OUT join("\t", ( ( map {$_ ."_sum"} @label) , ( map {$_ . "_avg"} @label ) ) ), "\n";

#for my $i_tmp (0..($file_num - 2)){
#	print OUT  join("\t", ( $label[$i_tmp] . "_sum",  $label[$i_tmp] . "_avg")) , "\t";
#}

#last one
#my $i_tmp = $file_num - 1;
#print OUT  join("\t", ( $label[$i_tmp] . "_sum",  $label[$i_tmp] . "_avg")) , "\n";

while (<IN>) {
	chomp;
	my @a = split "\t";
	my $chr = simple_chr ($a[0]);
	$chr = "Chr" . $chr;
	my ( $start, $end)  = @a[1..2];
	my $length = $end - $start + 1;
	my $interval = $chr . ":" . $start . "-" . $end ;
	my $cmd = "samtools depth -r $interval @bam_files_full";
	my @deps = `$cmd`;
	
	my @sums = ();
	
	foreach my $line (@deps){
		chomp $line;
		my @pts = split "\t", $line;
		
		for my $j (2..$#pts){
			$sums[$j - 2] += $pts[$j];
		}
	}
	
	my @avgs = map {sprintf("%.3f", $_ / $length)} @sums;
	
	if ( @deps == 0) {
		@sums = ("NA") x $file_num;
		@avgs = ("NA") x $file_num;
	}
	
	
	print OUT join("\t", (@a, ,@sums , @avgs )), "\n";	
}

close IN;
close(OUT);

exit;

#my @label = get_bam_label(@bam_files_simple);
sub get_bam_label{
	my (@tmp) = @_;
	my @lab = ();
	foreach my $item (@tmp){
		if ($item =~ /(\S+)\.bam$/) {
			push @lab, $1;
		}
		
	}
	return @lab;
}


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
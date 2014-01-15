#!/usr/bin/perl -w
# Copyright (C)  	2011-
# Program:			fastq_profile.pl
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		Mar 13, 2011
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	
# Description:		input clean_fastq file, count the length distribution 
#					and nucleitide distribution,the firs five and last five nucleitide
#					distribution
#**************************
# Version: 1.0
# 
#**************************
# e-mail:tangkai.ustc@gmail.com

my $version="1.0";
use strict;
use warnings;

use Getopt::Long;


my $debug = 0;
### Command line parameters -------------------------------------------------------------------
my $input	= "";	#input fastq file;
my $output	= "";	#output fasta file;
my %CMD;

print STDERR ("\n==================| $0 start |==========================================\n");
my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";

GetCom();

my ($total, $As, $Ts, $Cs, $Gs) = (0, 0, 0, 0, 0);
my %lens;
#my ($first_one, $first_two, $first_three, $first_four, $first_five) = (0, 0, 0, 0, 0);
#my ($last_one, $last_two, $last_three, $last_four, $last_five) = (0, 0, 0, 0, 0);
#my (%one, %two, %three, %four, %five);
my (%heads, %tails);
open(IN, $input) or die "cannot open $input:$!";
if (-e $output){die "$output exists!!!\n";}
open (OUT, ">$output") or die "cannot open $output:$!";
my $i = 0;

while(<IN>){
	$i++;
	if($i % 4 == 2){
		chomp;
		my $l = length($_);
		$lens{$l} ++;
		$total += $l;
		$As += ($_ =~ tr/Aa/Aa/);
		$Ts += ($_ =~ tr/Tt/Tt/);
		$Cs += ($_ =~ tr/Cc/Cc/);
		$Gs += ($_ =~ tr/Gg/Gg/);
		for (my $j = 0; $j <= 4; $j ++){
			my $head = substr($_, $j, 1);
			my $tail = substr($_, $l - 1 - $j, 1);
			$heads{$j}->{$head} ++;
			$tails{$j}->{$tail} ++;
		}
	}
}

close(IN);
my $max = 0;
my $min = 99999;
foreach my $l (sort {$a <=> $b} keys %lens){
	if($l > $max) {$max = $l}
	if($l < $min) {$min = $l}
}

print OUT join("\t", ("total_base", "A", "T", "C", "G")), "\n";
print OUT join("\t", ($total, $As, $Ts, $Cs, $Gs)), "\n\n";

print OUT "reads_length\t", "reads_number","\n";
for(my $j = $min; $j <= $max; $j++){
	my $out = 0;
	if(defined $lens{$j}){
		$out = $lens{$j}
	}
	print OUT $j,"\t", $out, "\n";
}
print OUT "\n";
print OUT join("\t", ("position","base", "head", "tail")), "\n";
for (my $j = 0; $j <= 4; $j++){
	foreach my $base ("A","T","C","G"){
		my $out_head = 0;
		my $out_tail = 0;
		if (defined $heads{$j}->{$base}){$out_head = $heads{$j}->{$base}}
		if (defined $tails{$j}->{$base}){$out_tail = $tails{$j}->{$base}}
		print OUT join("\t", ($j + 1, $base, $out_head, $out_tail)), "\n";
	}
}
close(OUT);


sub_end_program();
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}

############################################################################################################
######################                  Read command line parameters
############################################################################################################
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory:
--input		STRING
--output	STRING
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "input=s", "output=s");
	
	# Mandatory params
	die("Please specify input file\n") unless defined($CMD{input});
	die("Please specify output file\n") unless defined($CMD{output});
	$input	= $CMD{input};	#input fastq file;
	$output= $CMD{output};	#output fasta file;
}


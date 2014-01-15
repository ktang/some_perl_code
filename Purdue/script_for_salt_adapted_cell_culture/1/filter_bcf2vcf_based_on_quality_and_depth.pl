#!/usr/bin/perl -w
#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <in_bcf> <output_vcf>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

#open(IN, $input) or die "cannot open $input: $!";
#open( IN, "samtools view -h $input |") or die;
#open( OUT,  "| samtools view -b -S -  > $output") or die;

open (IN,  "bcftools view $input |") or die;
die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

#0	1	2	3	4	5	6	7	8	9
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ../MAPQ20_bam_ln/2/Col0_WT_0_labelAsrdd_0_bowtie2_sorted_MAPQ20.bam
#chr1    59214   .       C       A       7.8     .       DP=14;VDB=0.0216;AF1=0.5;AC1=1;DP4=2,2,0,4;MQ=41;FQ=10.4;PV4=0.43,8.8e-05,0,1   PL      37,0,113





while (<IN>) {
#	chomp;
	my @t = split;
	if (/^#/) {
		print OUT ;
		next;
	}
	next if ($t[4] eq '.'); # skip non-var sites
	next if ($t[0] =~ /chrc/i or $t[0] =~ /chrm/i or  $t[0] =~ /Pt/i  or  $t[0] =~ /Mt/i  );
	my $flt = 0;
	# parse annotations
	my ($dp, $mq, $dp_alt) = (-1, -1, -1);
	
	if ($t[7] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/i) {
	#  $dp = $1 + $2 + $3 + $4;
	  $dp_alt = $3 + $4;
	}
	if ($t[7] =~ /DP=(\d+)/i) {
	  $dp = $1;
	}
	my $qual = $t[5];
	
#1, only retain chr1-5;
#2, SBS (none INDEL)
#score >= 20 coverage 8-100;
#INDEL
# coverage 5-100
	# check if the site is a SNP
	if( $t[7] =~ /^INDEL/ ){
		print OUT if ( $dp >= 5 and $dp <= 100 );
		
	}elsif( $t[7] =~ /^DP=/  ){
		print OUT if (  $dp >= 8 and $dp <= 100 and $qual >= 20);
	}else{
		die $_;
	}

	

}

close(IN);
close(OUT);

exit;
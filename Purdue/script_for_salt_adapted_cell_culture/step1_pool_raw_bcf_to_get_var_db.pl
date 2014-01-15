#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

# input and output is std
# cat *vcf | perl $0 > output.txt

use strict;
use File::Spec;
# chr1    3       .       CTAAACCCTAAACCCTAAACCCTAAACC    CTAAACCCTAAACCCTAAACCCTAAACCCTAAACC     34.2    .       INDEL;DP=3;VDB=0.0060;AF1=1;AC1=2;DP4=0,0,2,0;MQ=40;FQ=-40.5    GT:PL:GQ        1/1:73,6,0:10
# chr1    59214   .       C       A       7.8     .       DP=14;VDB=0.0216;AF1=0.5;AC1=1;DP4=2,2,0,4;MQ=41;FQ=10.4;PV4=0.43,8.8e-05,0,1   GT:PL:GQ        0/1:37,0,113:39
# 0       1       2        3       4       5      6       7

my $usage = "$0 \n <indir_raw_bcf> <output>\n\n";
die $usage unless (@ARGV == 2);
my $indir = shift or die;
my $output = shift or die;

die unless ( -d $indir);
die if(-e $output);

open(OUT, ">>$output") or die; 
print OUT join("\t",("chr", "pos") ), "\n";
	       
opendir (DIR, $indir) or die ;
my @files = grep /raw\.bcf$/ , readdir DIR;
closedir DIR;

my %records;


foreach my $file(@files){
	my $input = File::Spec->catfile($indir, $file);
	die $file unless (-e $input);
	open(IN, "bcftools view $input |") or die;
	
	my $i = 0;
	
	while(<IN>){
		next if (/^#/);
		chomp;
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		if ( $chr =~ /\d/) {
			$i++;
			$chr =~ s/chr//i;
			$records{$chr}->{$pos} = 1;

		}
		
		
	}
	print STDERR join("\t", ($file, $i)), "\n\n";
	close IN;
}



foreach my $chr (sort keys %records){
	foreach my $pos(sort {$a <=> $b} keys %{$records{$chr}}){
		#die unless( defined  $refs{$chr}->{$pos});
		print OUT join("\t", ($chr, $pos )), "\n";
	}
}
close OUT;
exit;


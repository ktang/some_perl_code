#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

# input and output is std
# cat *vcf | perl $0 > output.txt

use strict;
# chr1    3       .       CTAAACCCTAAACCCTAAACCCTAAACC    CTAAACCCTAAACCCTAAACCCTAAACCCTAAACC     34.2    .       INDEL;DP=3;VDB=0.0060;AF1=1;AC1=2;DP4=0,0,2,0;MQ=40;FQ=-40.5    GT:PL:GQ        1/1:73,6,0:10
# chr1    59214   .       C       A       7.8     .       DP=14;VDB=0.0216;AF1=0.5;AC1=1;DP4=2,2,0,4;MQ=41;FQ=10.4;PV4=0.43,8.8e-05,0,1   GT:PL:GQ        0/1:37,0,113:39
# 0       1       2        3       4       5      6       7

my %records;

my %refs;

while (<>){
	chomp;
	my @a = split "\t";
	my $indel = 0;
	my ($chr, $pos) = @a[0..1];
	my $snp = $a[4];
	
	my $ref = substr ($a[3], 0, 1);
	
	if($a[7] =~ /^INDEL/){
		$indel = 1;
	}
	
	if (defined $refs{$chr}->{$pos}  ) {
		
		if ( $refs{$chr}->{$pos} ne $ref ){
			print STDERR $_, "\n";
			print STDERR $_, "\n";
			print STDERR $ref, "\t", $refs{$chr}->{$pos}, "\n";
			die;
		}
	}else{
		$refs{$chr}->{$pos} =  $ref;
	}
	
	
	if( !defined  $records{$chr}->{$pos} ){
		if( $indel){
			$records{$chr}->{$pos}  = "INDEL";
		}else{
			print STDERR $_, "\n" if (length($snp) != 1);
			$records{$chr}->{$pos}  = $snp;
		}
	}else{
		if( $records{$chr}->{$pos} =~ /INDEL/ and $indel ){
			next;
		}elsif( $records{$chr}->{$pos} =~ /INDEL/ and !$indel ){
			 $records{$chr}->{$pos}  =  $records{$chr}->{$pos}  . "," . $snp;
		}else{
			 $records{$chr}->{$pos}  =  $records{$chr}->{$pos}  . "," . $snp;
		}
		
	}
}

foreach my $chr (sort keys %records){
	foreach my $pos(sort {$a <=> $b} keys %{$records{$chr}}){
		die unless( defined  $refs{$chr}->{$pos});
		print join("\t", ($chr, $pos, $refs{$chr}->{$pos}, $records{$chr}->{$pos} )), "\n";
	}
}
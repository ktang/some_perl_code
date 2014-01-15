#!/usr/bin/perl -w

# annotate_individual_base_v2.0.pl
# April 18, 2013
# Purpose: On Apr18, for Xiaohong's Cell Culture project
# for each variation loci
# decide which feature is it
# 1 => CDS
# 2 => TE
# 3 => non_coding_RNA
# 4 => UTR
# 5 => intron
# 6 => pseudogene
# undef => IG

use strict;
use File::Spec;

my $debug = 0 ;
if($debug){
	print STDERR "debug:$debug\n";
}


# should be GFF file
my $gff_dir = "/Users/tang58/DataBase/TAIR10/GFF/Kai/annotate_individual_base";
my $CDS_file		= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_CDS.gff" ) ;
my $TE_file		= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_for_TE.gff" ) ;
my $non_coding_file	= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_for_non_coding_RNA.gff" ) ;
my $UTR_file		= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_for_UTR.gff" ) ;
my $intron_file		= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_intron_Kai.gff" ) ;
my $pseudogene_file	= File::Spec->catfile ($gff_dir , "TAIR10_GFF3_for_pseudogene.gff" ) ;
#my $intergenic_file	= File::Spec->catfile ($gff_dir , "" ) ;

die unless ( -e $CDS_file);
die unless ( -e $TE_file );
die unless ( -e $non_coding_file);
die unless ( -e $UTR_file);
die unless ( -e $intron_file);
die unless ( -e $pseudogene_file);

# input format
# chr pos

my $usage = "$0 \n <input> <output> \n\n";
die $usage  unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if (-e $output);

my @lists;
#my %records;

my %code_types = (1 => "CDS",
		  2 => "TE",
		  3 => "ncRNA",
		  4 => "UTR",
		  5 => "intron",
		  6 => "pseudogene"
	#	  7 => "IG"
		  );

read_input ($input, \@lists);

for my $chr (1..5) {
	print STDERR "chr" , $chr, "\n";
	my @pos_code = ( );
	my @pos_type = ( );
	
	#record( $chr, \@pos_code,  1..6 , exon;intron.., $file );
	# order is very important
# 1 => CDS
# 2 => TE
# 3 => non_coding_RNA
# 4 => UTR
# 5 => intron
# 6 => pseudogene
# 7 undef => IG

	record( $chr, \@pos_code, 6 ,  	$pseudogene_file);
	record( $chr, \@pos_code, 5 , 	$intron_file);
	record( $chr, \@pos_code, 4 , 	$UTR_file);
	record( $chr, \@pos_code, 3 , 	$non_coding_file);
	record( $chr, \@pos_code, 2 , 	$TE_file);
	record( $chr, \@pos_code, 1 , 	$CDS_file);
	
	put_type ( $chr, \@pos_code, \@lists, \%code_types );
	
}
output_file ($output, \@lists);

exit;

#output_file ($output, \@lists);
sub output_file{
	my ($file, $list_ref_a) = @_;
	open(OUT, ">>$file") or die;
	
	my $head = $list_ref_a ->[0];
	print OUT join("\t", ($head, "type", "type_code")), "\n";
	for my $i( 1..(  scalar(@{$list_ref_a}) -1  )   ){
		print OUT $list_ref_a ->[$i], "\n";
	}
	close OUT;
}

# put_type ( $chr, \@pos_code, \@lists, \%code_types );
# put last col as type and code;
# @lists has head
sub put_type{
	my (  $chr_num, $pos_code_ref_a,  $list_ref_a , $code_types_ref_h) = @_;
	my $last_index = scalar (@{$list_ref_a }) - 1;
	for my $i ( 1..$last_index){
		my $this = $list_ref_a->[$i];
		my @a = split "\t", $this;
		my ($chr, $pos ) = @a[0..1];
		$chr =~ s/chr//i;
		next unless ($chr eq $chr_num);
	#	my ($code, $type) = (7, "IG");
		my ( $type, $code) = ("IG", 7);
		
		if (defined $pos_code_ref_a->[$pos]){
			$code = $pos_code_ref_a->[$pos];
			if( $code =~ /;/){
				my @tmp = split ";", $code;
				$type = $code_types_ref_h->{$tmp[0]};
				for my $k (1..$#tmp){
					$type .=  ";". $code_types_ref_h->{$tmp[$k]};
				}
				
			}else{
				die $this unless (defined $code_types_ref_h->{$code}) ;
				$type = $code_types_ref_h->{$code};
			}
		}
		
		$list_ref_a->[$i] = join("\t", (@a, $type, $code));
		
	}
}

#GFF
#Chr1    TAIR9   CDS     3760    3913    .       +       0       Parent=AT1G01010.1,AT1G01010.1-Protein;
# record ( $chr, \@pos_code,  ,  );
sub record{
	my ( $chr_num, $pos_code_ref_a,  $code_sub,  $file ) = @_;
	die unless (-e $file);
	
	print STDERR $file, "\n";
	
	open(IN, $file) or die;
	while (<IN>){
		chomp;
		my @a = split "\t";
		my ($chr, $start, $end ) = @a[0,3,4];
		$chr =~ s/chr//i;
		next unless ( $chr eq  $chr_num);
		
		for my $i( $start..$end ){
			if (defined $pos_code_ref_a->[$i]){
				$pos_code_ref_a->[$i] .= ";". $code_sub;
			}else{
				$pos_code_ref_a->[$i] = $code_sub;
			}
		}
	}
	close IN;
}

#read_input ($input, \%lists, \%records);
sub read_input{
	my ($file, $list_ref_a) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		$list_ref_a->[$i] = $_;
		next if ($i == 0); # skip head
		my ($chr, $pos) = @a[0..1];
		$chr =~ s/chr//i;
	}
	close IN;
}
#!/usr/bin/perl -w

# miRNA length 20-24
# ReadMe
# for each mature RNA find its hairpin, include +- 5bp
# find all possible 20-24bp seq


#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;
use Bio::SeqIO;

my $bp = 5;

my $dir = "/Users/tang58/DataBase/miRBase/Release_19";
die unless ( -d $dir);
my $ath_mature_fa = File::Spec->catfile ( $dir, "ath_mature_r19_T.fa");
my $ath_hairpin_fa = File::Spec->catfile ( $dir, "ath_hairpin_r19_T_oneline.fa" );
die unless ( -e $ath_mature_fa);
die unless ( -e $ath_hairpin_fa);

my $id_dir = "/Users/tang58/DataBase/miRBase/Release_19/Kai";

my $dot_file = File::Spec->catfile ( $id_dir, "ath_id_with_dot.txt");
my $p3_file  = File::Spec->catfile ( $id_dir, "ath_id_with_3p.txt");
my $p5_file   = File::Spec->catfile ( $id_dir, "ath_id_with_5p.txt");

die unless ( -e $dot_file);
die unless ( -e $p3_file );
die unless ( -e $p5_file );

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <input> <output>\n\n";
#die $usage unless(@ARGV == 2);

my $usage = "$0 \n <output>\n\n";
die $usage unless(@ARGV == 1);

my $output = shift or die;

die if (-e $output);

if($debug) {
	exit;
}

my $seqin = Bio::SeqIO->new(-file=>$ath_mature_fa, -format=>'fasta');
my %miRNAs;
my %miRNA_ids;
while(my $seq = $seqin->next_seq){
	my $seqstr = $seq->seq;
#	$seqstr =~ tr/U/T/;
	my $id =  lc( $seq->id);
	
	my $simple_id = "";
	if( $id =~ /(ath-mir\d+[a-z]?)[\.?|-]?/){
		$simple_id = $1;
		die $simple_id if( $simple_id =~ /p/ or $simple_id =~ /\./ );
	}else{
		die "simple id: ", $id;
	}
	
	
	$miRNA_ids{$id} = [$simple_id, $seqstr];
	if(!defined $miRNAs{$seqstr}){
		$miRNAs{$seqstr} = $id;
	}else{
		$miRNAs{$seqstr} = $miRNAs{$seqstr} . ";" . $id;
	}
}
$seqin->close;

my %extend_seq;
my %hairpin_seq;

read_hairpin_file( $ath_hairpin_fa , \%hairpin_seq );

extract_all_extend_seq( \%hairpin_seq,  \%extend_seq, \%miRNAs, \%miRNA_ids);



open (OUT, ">>$output") or die;

foreach my $key (sort keys %extend_seq ){
	print  OUT join("\t",  ( $key, $extend_seq{$key} )), "\n";
}

close OUT;
exit;
sub read_hairpin_file{
	my ($fa_file, $ref) = @_;
	die unless ( -e $fa_file);
	my $seqin = Bio::SeqIO->new(-file=>$fa_file, -format=>'fasta');
	while(my $seq = $seqin->next_seq){
		my $seqstr = $seq->seq;
		my $id =  lc( $seq->id);
		$ref->{$id} = $seqstr;
	}
	$seqin->close;
}

# extract_all_extend_seq( \%hairpin_seq,  \%extend_seq, \%miRNAs, \%miRNA_ids);
sub extract_all_extend_seq{
	my ( $hairpin_seq_ref,  $extend_seq_ref, $miRNAs_ref, $miRNA_ids_ref ) = @_;
	foreach my $raw_id (sort keys %{$miRNA_ids_ref}){
		my ($simple_id, $mature_seq ) = @{ $miRNA_ids_ref->{$raw_id} };
		die $simple_id unless (defined $hairpin_seq_ref->{$simple_id} );
		
		my $seqstr = $hairpin_seq_ref->{$simple_id};
		
		my $length_miRNA = length ( $mature_seq );
		my $length_hairping = length ( $seqstr ) ;
		
		if( $seqstr =~ /$mature_seq/ ){
			my $match_index = $-[0];
			my $time = ( $seqstr =~ s/$mature_seq/$mature_seq/g );
			if ($time > 1){
				die "appear more than one ", $raw_id;
			}
			my $start_index =  ( $match_index - $bp >= 0  ) ? (  $match_index - $bp ) : ( 0 ) ;
			my $end_index   =  ( $match_index + $length_miRNA - 1 + $bp <  $length_hairping - 1 ) ?
					   ( $match_index + $length_miRNA - 1 + $bp) : ($length_hairping - 1 );
		
			my $extend_seq_len = $end_index - $start_index + 1;
			my $extend_seq = substr( $seqstr ,  $start_index ,  $extend_seq_len);
			die unless ( $extend_seq =~ /$mature_seq/);
		
			for my $i (20..24){
				for my $j (0..( $extend_seq_len - $i)){
					my $tmp = substr( $extend_seq, $j, $i );
					next if ( defined $miRNAs_ref->{$tmp} );
					if( !defined $extend_seq_ref->{$tmp}){
						$extend_seq_ref->{$tmp} = $simple_id;
					}else{
						$extend_seq_ref->{$tmp} =  $extend_seq_ref->{$tmp}  . ";" . $simple_id;
					}
				}
			}
		}else{
			die "no match",  $raw_id;
		}
	}
}
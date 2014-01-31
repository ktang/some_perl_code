#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use strict;
use File::Spec;

my $debug = 0;

my %chr_len = ("1"=>30427671, "2"=>19698289, "3"=>23459830, "4"=>18585056, "5"=>26975502);

my $win_size = 1000;
my $sliding_bp = 500;

my $strand = "+";

my %num_win;
my %cumsum_win_SE;


if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <output_pre> \n\n";
die $usage unless(@ARGV == 1);

#my $input = shift or die;
my $outpre = shift or die;

#die unless (-e $input);
my $output = $outpre . "_WinSize" . $win_size. "bp_sliding" . $sliding_bp . "bp_SAF.txt";
print STDERR $output, "\n\n";
die if( -e $output);

#open(IN, $input) or die "cannot open $input: $!";



Kai_Module::cal_win_num( \%num_win, \%chr_len, $win_size, $sliding_bp);
Kai_Module::hash_cumsum_SE ( \%cumsum_win_SE,  \%num_win );


#1	1	60856	60856
#2	60857	100253	39397
#3	100254	147173	46920
#4	147174	184344	37171
#5	184345	238296	53952


if ($debug) {
	foreach my $chr( sort keys %cumsum_win_SE){
		my ( $s_ind ,$e_ind ) = @{$cumsum_win_SE{$chr}};
		print STDERR join("\t", ($chr, $s_ind ,$e_ind, $e_ind - $s_ind +1)), "\n";
	}
	exit;#code
}

die if(-e $output);
if (!$debug) {
	open(OUT, ">$output") or die "cannot open $output: $!";
#code
}

foreach my $chr( sort keys %cumsum_win_SE){
	my ( $s_ind ,$e_ind ) = @{$cumsum_win_SE{$chr}};
	#print STDERR join("\t", ($chr, $s_ind ,$e_ind, $e_ind - $s_ind +1)), "\n";
	for my $id ( $s_ind.. ( $e_ind - 1) ){
		my $start = ($id - $s_ind) * $sliding_bp + 1;
		my $end = $start + $win_size - 1;
		print  OUT join("\t", ($id, $chr, $start, $end , $strand)), "\n"
	}
	my $id = $e_ind;
	my $start = ($id - $s_ind) * $sliding_bp + 1;
	my $end = $chr_len{$chr};
	print  OUT join("\t", ($id, $chr, $start, $end , $strand)), "\n";
}
#close(IN);
close(OUT);

exit;


#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);


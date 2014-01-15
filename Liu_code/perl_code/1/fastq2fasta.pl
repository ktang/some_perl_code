#!/usr/bin/perl -w
# take illumina data in fastq format and output sequence data in fasta format
# without considering quality
use strict;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "$0 <fastq file> <output seq file> [<output qual file>]";
die $usage unless (@ARGV >=2);
my ($input, $output) = @ARGV[0..1];
my $qual_file;
if(@ARGV >=3){
	$qual_file = $ARGV[2];
	
}
my $out = Bio::SeqIO->new(-file=>">$output", -format=>"fasta");
open(IN, $input) or die "Cannot open $input:$!";
my $qual_out;
if(defined $qual_file){
	open $qual_out, ">$qual_file" or die "cannot open $qual_file:$!";
}
my $line = 0;
my ($seq_name, $seq_str);
my ($num_yes, $num_no) = (0,0);
while(<IN>){
	$line++;
	chomp;
	my $n = $line % 4;
	if($n == 1){
		if(/\@(\d+)\:(\d+)\:(\d+)\:(\d+)\:([Y|N])/){
			$seq_name = join("_", ($1, $2, $3, $4, $5));
		}else{
			die "On line $line, seq name $_ does not match pattern\n";
		}
		if($5 eq 'Y'){
			$num_yes++;
		}
		if($5 eq 'N'){
			$num_no++;
		}
	}
	if($n == 2){
		$seq_str = $_;
		my $seq = Bio::Seq->new(-seq => $seq_str, -id=>$seq_name);
		$out->write_seq($seq);
	}
	if((defined $qual_file) && ($n == 0)){
		my @qual_vals;
		foreach my $i(1..length($_)){
			my $ch = substr $_, ($i - 1), 1;
			#my $q = int((10 * log(1+10**(ord($ch)-64)/10.0))/log(10));
			my $q = ord($ch) - 33;
			push @qual_vals, $q;
		}
		print $qual_out ">$seq_name\n";
		print $qual_out join(" ", @qual_vals), "\n";
	}
}
$out->close;
if(defined $qual_file){
	$qual_out->close;
}
print STDERR "Number of 'Y' seqs: $num_yes, Number of 'N' seqs: $num_no\n";	

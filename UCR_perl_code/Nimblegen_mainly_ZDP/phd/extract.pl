#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Bio::DB::Fasta;

# This program will count the number of sequenes in the files
# and count the number of bases in the whole file and count the 
# total numbers of A,C,G,T

#my $inputfile = shift @ARGV;

#phd_18_15_ave_minus_C2_L2L3_merged_positive_trimmean.bed
my $usage = "$0 <inFasta> <pos> <outFile>";

die $usage unless (@ARGV == 3);

my ($inputfile, $pos_file, $output) = @ARGV[0..2];

my $db = Bio::DB::Fasta->new($inputfile);

open (POS,$pos_file) 
	or die "cannot open $pos_file:$!";


#my @locations = (['scer_S288C_chrVII',100,150],
 #                ['scer_S288C_chrVII',1000,1050]);

my $out = Bio::SeqIO->new( -format=>'fasta',
			-file => ">$output");



my @locations ;

while (<POS>){
	chomp;
	my @pts = split "\t";
	$pts[0] =~ tr/c/C/;
	push @locations,[$pts[0],$pts[1],$pts[2]];
}
	

for my $loc (@locations)
{
	my($seqid,$start,$stop)= @$loc;
	my $seq = $db->seq($seqid,$start=>$stop);
	my $seqobj = Bio::Seq-> new(-seq => $seq,
			   # -desc => "$seqid:$start..$stop",
		        -id =>"$seqid".":"."$start"."_"."$stop");
	$out->write_seq($seqobj);
}

=head2
print "input file is $inputfile\n";

#my $filehandle;
#open($filehandle=>$inputfile);
#my next_data = <$filehandle>

my $in = Bio::SeqIO-> new(-file => $inputfile,
			  -format=> "fasta");
my $seq = $in -> next_seq;
#print $seq->

my $num_sequences;#we'll count the No. of seqs in the file here
my $num_bases;# we'll count total lentth of the database;
my %bases;

while ( $seq = $in-> next_seq)
{
  $num_sequences++;
  $num_bases += $seq->length;
  my $str = $seq->seq; #get a string
 
$base{"A"} += ($str =~ tr/Aa/Aa/);

#  my $a_count = ($str =~ tr/Aa/Aa/);
#  print $a_count,"\n";
 # last;
}

#printf "A = %.2f C = %.2f 
=cut
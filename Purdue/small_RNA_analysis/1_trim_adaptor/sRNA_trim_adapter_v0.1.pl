#!/usr/bin/perl -w
use strict;
use File::Spec;

#v0.1
# as I found some fq file have barcode seq at the first Xbp
# this version can specify the trim first N bp

#sRNA_trim_adapter.pl
# trim adapter for sRNA-seq
# input format can be fa fa_gz fq fq_gz

#output is fa file

my $debug = 0;

my $min_length = 18;
my $max_length = 31;

#fa 0
#fq 1
#gz 1
# non-gz 0
my %input_formats = (	fa	=> [0, 0],
			fa_gz 	=> [0, 1],
			fq 	=> [1, 0],
			fq_gz 	=> [1, 1]
		);

#my $usage = "$0 \n <input_coordinate_data> <db_for_extraction> <label> <output>  \n\n";
#die $usage unless(@ARGV == 4);

#AGATC
#GGAAG
#AGCAC
#ACGTC

#AGATCGGAAGAGCACACGTC

#my $usage = "$0 \n <adapter_seq> <input_format> <input_file> <outdir> <outpre>  \n\n";
#die $usage unless(@ARGV == 5);

my $usage = "$0 \n <adapter_seq> <input_format> <input_file> <outdir> <outpre> <trim5_bp> \n\n";
die $usage unless(@ARGV == 6);


my $adapter_seq = shift or die;
my $input_format = shift or die;
my $input 	= shift or die;

my $outdir = shift or die;
my $outpre = shift or die;

my $trim5_bp = shift ;
die unless ( $trim5_bp >= 0);

die unless  ( -e $input);
die unless (-d $outdir);

unless (defined $input_formats{$input_format}){
	print STDERR "format should be \n ";
	print STDERR (keys %input_formats);
	die;
}
 

my $out_fa = File::Spec->catfile($outdir, $outpre . "_clean_" .$min_length . "_" . $max_length . "bp.fa.gz" );
my $out_stat = File::Spec->catfile($outdir, $outpre . "_trim_stat.txt" );

die if ( -e $out_fa ) ; 
die if ( -e $out_stat) ;

print STDERR "ada = $adapter_seq\n\n";

if ($debug) {
	exit;#code
}

my ($fq_or_fa, $gz_or_not) = @{$input_formats{$input_format}};

read_fx($fq_or_fa, $gz_or_not ,$input, $out_fa, $out_stat, $adapter_seq, $min_length, $max_length, $trim5_bp);

=head
if ( $input_format eq "fa" ) {
	read_fa(0 ,$input, $out_fa, $out_stat, $adapter_seq, $min_length, $max_length, $trim5_bp);
}elsif($input_format eq "fa_gz"){
	read_fa(1 ,$input, $out_fa, $out_stat, $adapter_seq, $min_length, $max_length, $trim5_bp);
}elsif($input_format eq "fq"){
	read_fq(0 ,$input, $out_fa, $out_stat, $adapter_seq, $min_length, $max_length, $trim5_bp);
}elsif( $input_format eq "fq_gz" ){
	read_fq(1, $input, $out_fa, $out_stat, $adapter_seq, $min_length, $max_length, $trim5_bp)
}else{
	die ;
}
=cut

exit;

#read_fa
#read_fq
sub read_fx{
	my ($fq_label, $gz_label, $input_sub, $out_fa_sub, $out_stat_sub, $adapter_seq_sub, $min_length_sub, $max_length_sub, $trim5_bp_sub  ) = @_;
	die unless ( -e $input_sub);
	die if (-e $out_fa_sub);
	die if (-e $out_stat_sub);
	if ( $gz_label ) {
		die "file format" unless ( $input_sub =~ /gz$/ );
		open(IN, "zcat $input_sub |") or die;
	}else{
		die "file format" if ( $input_sub =~ /gz$/ );
		open(IN, $input_sub) or die;
	}
	open(OUT, "| gzip >$out_fa_sub") or die;
	open(STAT, ">$out_stat_sub") or die;
	
	my $total = 0;
	my $N_no_ada = 0;
	my $N_trimmed_but_N = 0;

	my %N_too_short;
	my %N_too_long;
	
	my %nums_accepted;
	while ( my $head = <IN>) {
		$total++;
		my $seq = <IN>;
		if ($fq_label) {
			my $plus = <IN>;
			my $qual =<IN>;#code
		}
		
		
		chomp $head;
		chomp $seq;
		
		$seq = substr($seq, $trim5_bp_sub);
		
		#if ( $seq =~ /^([ACTG]+)$adapter_seq_sub/) {
		if ( $seq =~ /(\S+)$adapter_seq_sub/) {
			my $trimmed_seq = $1;
			
		#	$trimmed_seq = substr($trimmed_seq, $trim5_bp_sub);
			
			my $length = length( $trimmed_seq );
			
			
			if ( $trimmed_seq =~ /N/) {
				$N_trimmed_but_N++;
				next;#code
			}elsif(	$length < $min_length_sub){
				#$N_too_short++;
				$N_too_short{$length}++;
				$N_too_short{sum}++;
				next;
			}elsif( $length >  $max_length_sub ){
				#$N_too_long++;
				$N_too_long{$length}++;
				$N_too_long{sum}++;
				next;
			}else{
				$nums_accepted{$length}++;
				#$head =~ s/@/>/;
				if ($head =~ /(\d+:\d+:\d+:\d+)\s*/) {
					my $out_head = ">" . $1;
					print OUT join("\n", ($out_head, $trimmed_seq)), "\n";
					next;
				}else{
					die  $head;
				}
				
			}
			
		}else{
			$N_no_ada++;
		}
	}
	close OUT;
	close IN;
	
	print STAT join("\t", ("type", "number")), "\n";
	
	my $sum = 0;
	for my $length ($min_length_sub..$max_length_sub ){
		my $n = 0;
		if (defined $nums_accepted{$length} ) {
			$n = $nums_accepted{$length};
		}
		$sum+=$n;
		print STAT join("\t", ($length, $n) ), "\n";
		
	}
	
	print STAT join("\t", ("sum", $sum) ), "\n\n\n";
	print STAT join("\t", ("Total",  $total) ), "\n";
	print STAT join("\t", ("no_adapter", $N_no_ada ) ), "\n";
	print STAT join("\t", ("trimmed_but_has_N", $N_trimmed_but_N) ), "\n";
	
	if (defined $N_too_short{sum} ) {
		print STAT join("\t", ("too_short", $N_too_short{sum} ) ), "\n";
		delete $N_too_short{sum};
		foreach my $i (sort {$a<=>$b} keys %N_too_short){
			print STAT join("\t", ($i, $N_too_short{$i} ) ), "\n";
		}
	}
	
	if (defined $N_too_long{sum} ) {
		print STAT join("\t", ("too_long", $N_too_long{sum} ) ), "\n";
		delete $N_too_long{sum};
		foreach my $i (sort {$a<=>$b} keys %N_too_long){
			print STAT join("\t", ($i, $N_too_long{$i} ) ), "\n";
		}
	}
	close STAT;
}

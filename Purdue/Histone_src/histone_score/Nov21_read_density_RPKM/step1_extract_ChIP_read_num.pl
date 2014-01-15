#!/usr/bin/perl -w

#Nov 21, 2013: v0.0
#step1_extract_ChIP_read_num.pl
# input DMR_bed_like list, Notice that the chr name of bam file(for color space) is
# Chr1-5, chloroplast, mitochondria
use strict;
use File::Spec;

my $debug = 0;
my $print_debug = 0;

if ( !$debug ) {
	$print_debug = 0;
}
if($debug){
	print STDERR "debug = 1\n\n";
}

#my $half_interval_length = 2000;
# order for bins: -25,-24,...,-1,0, 1,...25;
my $half_flanking_bin_num = 25; #100; #40;
my $bin_size = 200;#50;

my %types = ("bed"=>1, "coor" => 1);

my $usage = "$0 \n <input_list_DMR> <type: bed or coor> <bam_file>  <bam_label> <read_length> <outdir> <output_pre>\n\n";
die $usage unless(@ARGV == 7);

my $input	= shift or die;
my $type	= shift or die;
die "type: bed or coor" unless (defined $types{$type});
#my $bam_dir	= shift or die;
my $bam_file 	= shift or die;
my $bam_label	= shift or die;
my $read_length = shift or die;
my $outdir	= shift or die;
my $outpre	= shift or die;

#my $output	= shift or die;
#my $output = File::Spec->catfile($outdir, $outpre . "_half" . $half_interval_length . "bp" . "_bin" . $bin_size ."bp.txt");
my $output = File::Spec->catfile($outdir, $outpre . "_in_" . $bam_label . "_halfBinN" . $half_flanking_bin_num  . "_bin" . $bin_size ."bp_ReadsN.txt");

die unless (-e $input);
die if( -e $output);
#die unless (-d $bam_dir);

###########
# 1, read_list_and_get_interval
##############

#Chr1	3	2
open(IN, $input) or die "cannot open input";
open(OUT, ">>$output") or die"cannot output";
print OUT join("\t", ("loci", "sum", (-1 * $half_flanking_bin_num)..$half_flanking_bin_num )), "\n";
while (<IN>) {
	my @reads_num = ();
	chomp;
	my @a = split "\t";
	next if ($a[0] =~ /(^coor|^#?chr$)/i );
	my ($chr, $start, $end);
	my $loci;
	if ($type eq "bed") {
		($chr, $start, $end) = @a[0..2];
		$chr = simple_chr($chr);
		$loci = "$chr:$start-$end";
		$chr = add_chr($chr);
	}elsif($type eq "coor"){
		if ($a[0] =~ /(\S+):(\d+)-(\d+)/) {
			($chr, $start, $end) =  ($1, $2, $3);
			$chr = simple_chr($chr);
			$loci = "$chr:$start-$end";
			$chr = add_chr($chr);
		}
		
	}else{
		die "wrong type\n\n";
	}
	
	my $mid_point = int ( ($start + $end)/2);
	#print STDERR "mid_point: $mid_point\n";
	
	my $first_start = $mid_point - int( ($bin_size - 1) /2 ) - $half_flanking_bin_num * $bin_size + $bin_size - 1 - int ( ($read_length)/2);
	if($first_start <=0){
		print STDERR "$loci is close to start \n";
		next;
	}
#	----------|----------
#	         + ++
	for my $bin_order ( (-1 * $half_flanking_bin_num)..$half_flanking_bin_num ){
		my $interval_start ;
		if ($bin_order <=0 ) {
			$interval_start = $mid_point - int( ($bin_size - 1) /2 ) +  $bin_order * $bin_size ;# original start
		}else{
			 $interval_start = $mid_point + int( ($bin_size + 1) /2 ) + 1 + ( $bin_order - 1) * $bin_size ;
		}
		
		
		my $old_start = $interval_start;
		
		my $interval_end   = $interval_start + $bin_size - 1 - int ( ($read_length)/2); # exclude half read length
		if ($interval_start <= 0) {
			if ($interval_end > 0) {
				$interval_start = 1;
			}else{
				push @reads_num, 0;
			}			
		}else{
			$interval_start = $interval_start + int ( ($read_length+1)/2) - 1; # real_start to exclude half read length
			
		}
		my $interval = $chr . ":" . $interval_start . "-" . $interval_end ;
		my $cmd = "samtools view -c $bam_file $interval";
		if ($print_debug) {
			print STDERR $old_start, "\t",$cmd, "\n";
		}
		my $num = `$cmd`;
		chomp $num;
		push @reads_num, $num;		
	}
	if (@reads_num != 1+ 2* $half_flanking_bin_num){
		print STDERR join("\n", @reads_num), "\n";
		die;
	}
	
	my $sum = get_sum(@reads_num);
	
	print OUT join("\t", ( $loci, $sum, @reads_num ) ), "\n";
	
	#print STDERR "\n\n";
}

close IN;
close OUT;

exit;
sub get_sum{
	my @tmp = @_;
	my $s = 0;
	foreach my $i(@tmp){
		$s+=$i;
	}
	return $s;
}

# read_list ($input, \@DMRs_list,  $bin_size, $half_flanking_bin_num);
sub read_list{
	my ($file, $ref_a, $bin_size_sub, $half_flanking_bin_num_sub) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $h = <IN>;
	while ( <IN> ) {
		chomp;
		my @a = split "\t";
		my $chr = simple_chr ($a[0]);
		my ($start, $end ) = @a[1..2];
		my $mid_point = int( ($start + $end)/2 );
		my $interval_label = $chr . ":" . $start . "-" . $end ;
		my $total_length = $bin_size_sub * ( 2 * $half_flanking_bin_num_sub + 1 );
		my $extend_start = $mid_point - ( $half_flanking_bin_num_sub * $bin_size_sub + int( ($bin_size_sub - 1) /2 ) ) ;
		my $extend_end   =  $mid_point + ( $half_flanking_bin_num_sub * $bin_size_sub + int( ($bin_size_sub) /2 ) ) ;
		if ($extend_start <=0 ) {
			next;
		}
		my $interval_extend = "Chr" . $chr . ":" . $extend_start . "-" . $extend_end ;
		push @{$ref_a}, [$interval_label, $interval_extend ]
	}
	
	close IN;
}


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub add_chr{
	my ($chr) = @_;
	if ($chr =~ /^\d+$/) {
		return "Chr" . $chr;
	}elsif($chr eq "Mt" ){
		return "mitochondria";
	}elsif( $chr eq "Pt"){
		return "chloroplast"
	}else{
		die "wrong chr from add_chr: $chr\n\n";
	}
	
	#chloroplast, mitochondria
}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
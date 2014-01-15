#!/usr/bin/perl -w

#v0.1
# also for 18-32bp

################
#	step2.1: generate non-overlaping interval, HNA value (Lee 2012)
################################
#reads_for_normalization:XXX
#chr:s-e	hits_number	reads_num	HNA

use strict;
use File::Spec;
my $debug = 0;

my $TM = 10000000;

my %chr_len = (
	       "1"=>30427671, "2"=>19698289, "3"=>23459830, "4"=>18585056, "5"=>26975502,
	       "Mt" => 366924, "Pt" => 154478
	       );

#my $bin_size = 500;
my @bp_interval = 18..32;

# index start from 0

if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <input_bam> <label> <total_reads>  <output>\n\n";
#die $usage unless(@ARGV == 4);

# bed list format
#chr	start	end
my $usage = "$0 \n <input_bed_list> <input_bam> <label> <total_reads>  <output>\n\n";
die $usage unless(@ARGV == 5);


my $bed_file	= shift or die;
my $input	= shift or die;
my $label	= shift or die;
my $total_reads = shift or die;
my $output	= shift or die;

die unless ( -e $bed_file);
die unless ( -e $input);
die if     ( -e $output);

my @bed_list;
my %pos;
my $last_one_index = read_bed_file($bed_file, \@bed_list);
# head is index 0 and "\n" is chompped.
print STDERR "last: line_num: ", $last_one_index, "\n\n";
record_region_pos(\@bed_list, \%pos);


#open( IN, "samtools view -h $input |") or die;
open( IN, "samtools view $input |") or die;

die if(-e $output);
open(OUT, ">>$output") or die "cannot open $output: $!";

#chr:s-e	hits_number	reads_num	HNA

my $factor = $TM / $total_reads;

print OUT "#reads_for_normalization:", $total_reads, "\t",  "normalization_factor:", $factor, "\n";
print OUT join("\t", ( $bed_list[0], "reads_". $label , "hits_sum_decimal_". $label , "HNA_" . $label ) ), "\t";

for my $i ( $bp_interval[0] .. $bp_interval[ $#bp_interval - 1 ] ){
	print OUT join( "\t", ( $i . "bp_reads_". $label ,  $i . "bp_HNA_" . $label   )) , "\t";
}
print OUT join( "\t", ($bp_interval[ $#bp_interval] . "bp_reads_". $label ,  $bp_interval[ $#bp_interval] . "bp_HNA_" . $label  ) ) , "\n";
if($debug){
	print STDERR join( "\n", ( $input, $total_reads, $output)), "\n\n";
	exit;
}

my (%raw_reads_sum_integer, %hits_sum_decimal);
my (%raw_reads_sum_integer_len, %hits_sum_decimal_len);

while(<IN>){
	chomp;
	my @a = split "\t";
	
	my($chr, $start ) = ($a[2], $a[3]);
	die $_ unless (defined $chr_len{$chr});
	my $length = length ( $a[9] );
	my $end = $start + $length - 1;
	
	next unless (defined $pos{$chr}->[$start] or defined $pos{$chr}->[$end] );
		
	my $hit;
	if( $a[-1] =~ /NH:i:(\d+)$/ ){
		$hit = 1/$1;
	}else{
		die $_;
	}
	
	my $index = -1;
	
	if( defined $pos{$chr}->[$start]  ){
		$index =  $pos{$chr}->[$start];
	}elsif( defined $pos{$chr}->[$end]	){
		$index = $pos{$chr}->[$end];
	}else{
		#die $_, $start_index, $end_index;
		die $_;
	}
	$raw_reads_sum_integer{$chr}->[$index] += 1;
	$hits_sum_decimal{$chr}->[$index] += $hit;
	$raw_reads_sum_integer_len{$chr}->[$index]->[$length] += 1;
	$hits_sum_decimal_len{$chr}->[$index]->[$length] += $hit;
}
close(IN);

#$last_one_index
foreach my $i (1..$last_one_index){
	
	my @a = split "\t", $bed_list[ $i ];
	
	my ( $chr, $start, $end ) = @a[0..2];
	
	$chr = simple_chr($chr);
	$a[0] = $chr;
	
	my ($reads, $hits, $HNA) = ( 0, 0, 0 );
	if(defined  $raw_reads_sum_integer{$chr}->[$i]){ $reads = $raw_reads_sum_integer{$chr}->[$i] 	}
	if(defined $hits_sum_decimal{$chr}->[$i] ) { $hits =  $hits_sum_decimal{$chr}->[$i]}
	$HNA = $hits * $factor;
	
	print OUT join("\t", (@a, $reads, $hits, $HNA)), "\t";

	my $total = 0;
		
		#for my $len ( $bp_interval[0] .. $bp_interval[ $#bp_interval - 1 ] ){
		for my $len ( $bp_interval[0] .. $bp_interval[ $#bp_interval - 1 ] ){
			my ( $reads_sub, $hits_sub, $HNA_sub ) = (0, 0);
			if(defined  $raw_reads_sum_integer_len{$chr}->[$i]->[$len]  ){ $reads_sub = $raw_reads_sum_integer_len{$chr}->[$i]->[$len] 	}
			if(defined  $hits_sum_decimal_len{$chr}->[$i]->[$len]       ){ $hits_sub  = $hits_sum_decimal_len{$chr}->[$i]->[$len]       }
			$HNA_sub = $hits_sub * $factor;
			$total += $reads_sub;
			print OUT join("\t", ( $reads_sub, $HNA_sub )), "\t";
		}
		my $len = $bp_interval[ $#bp_interval];
		my ( $reads_sub, $hits_sub, $HNA_sub ) = (0, 0);
		if(defined  $raw_reads_sum_integer_len{$chr}->[$i]->[$len]  ){ $reads_sub = $raw_reads_sum_integer_len{$chr}->[$i]->[$len] 	}
		if(defined  $hits_sum_decimal_len{$chr}->[$i]->[$len]       ){ $hits_sub  = $hits_sum_decimal_len{$chr}->[$i]->[$len]       }
		$HNA_sub = $hits_sub * $factor;
		print OUT join("\t", ( $reads_sub, $HNA_sub )), "\n";
		$total += $reads_sub;
		
	die $chr, $start, $end  if ($total != $reads);
}



close(OUT);

exit;

sub read_bed_file{
	my ($file, $ref) = @_;
	die unless (-e $file);
	
	open(IN, $file) or die;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		$ref->[$i] = $_;
	}
	close(IN);
	
	return $i;
}


sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ( $chr, $start, $end ) = @a[0..2];
		
		$chr = simple_chr($chr);
		
		die join("\t", ( $chr, $start, $end) ) unless (defined $chr_len{$chr});

		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i; #record index
		}
	}
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
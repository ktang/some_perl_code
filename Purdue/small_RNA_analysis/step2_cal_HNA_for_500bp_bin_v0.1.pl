#!/usr/bin/perl -w

#v0.1
# also for 18-32bp

################
#	step2: generate table of 500bp-bin, HNA value (Lee 2012)
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

my $bin_size = 500;
my @bp_interval = 18..32;

# index start from 0

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <label> <total_reads>  <output (NoOver_strc_RNA_500bp_HNA_all/uniq.txt)>\n\n";
die $usage unless(@ARGV == 4);

my $input = shift or die;
my $label = shift or die;
my $total_reads = shift or die;
my $output = shift or die;

die unless (-e $input);
die if     ( -e $output);



#open( IN, "samtools view -h $input |") or die;
open( IN, "samtools view $input |") or die;

die if(-e $output);
open(OUT, ">>$output") or die "cannot open $output: $!";

#chr:s-e	hits_number	reads_num	HNA

my $factor = $TM / $total_reads;

print OUT "#reads_for_normalization:", $total_reads, "\t",  "normalization_factor:", $factor, "\n";

#HWI-ST531R:169:D1549ACXX:7:1101:3471:1999       16      2       16539980        255     21M     *       0       0       AAGTATCATCATTCGCTTGGA   IIIIIIIIIIIIIIIIIIIII   XA:i:0  MD:Z:21 NM:i:0  NH:i:1
#HWI-ST531R:169:D1549ACXX:7:1101:2827:1997       16      3       12098843        255     21M     *       0       0       AGGCCTCTACGAATTCATGAT   IIIIIIIIIIIIIIIIIIIII   XA:i:0  MD:Z:21 NM:i:0  NH:i:1
#HWI-ST531R:169:D1549ACXX:7:1101:5870:1984       0       2       8658085 	255     23M     *       0       0       ACGGAATAATGTAAAACTGTACA IIIIIIIIIIIIIIIIIIIIIII XA:i:0  MD:Z:23 NM:i:0  NH:i:1
#			0			 1	 2	  3	 	 4	5	6	7	8	 9				10	
#print OUT join("\t", ( "coordinate", "raw_reads_sum_integer". $label , "hits_sum_decimal". $label , "HNA" . $label ) ), "\n";

print OUT join("\t", ( "coordinate", "reads_". $label , "hits_sum_decimal_". $label , "HNA_" . $label ) ), "\t";
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
	
	my $hit;
	if( $a[-1] =~ /NH:i:(\d+)$/ ){
		$hit = 1/$1;
	}else{
		die $_;
	}
		
	my $start_index =  int (($start - 1) / $bin_size );
	my $end_index = int (( $end - 1) / $bin_size );
	if($start_index == $end_index){
		$raw_reads_sum_integer{$chr}->[$start_index] += 1;
		$hits_sum_decimal{$chr}->[$start_index] += $hit;
		
		$raw_reads_sum_integer_len{$chr}->[$start_index]->[$length] += 1;
		$hits_sum_decimal_len{$chr}->[$start_index]->[$length] += $hit;
		
	}elsif( $start_index ==  $end_index - 1){
	
		$raw_reads_sum_integer{$chr}->[$start_index] += 1;
		$hits_sum_decimal{$chr}->[$start_index] += $hit;
		$raw_reads_sum_integer{$chr}->[$end_index] += 1;
		$hits_sum_decimal{$chr}->[$end_index] += $hit;
		
		$raw_reads_sum_integer_len{$chr}->[$start_index]->[$length] += 1;
		$hits_sum_decimal_len{$chr}->[$start_index]->[$length] += $hit;
		$raw_reads_sum_integer_len{$chr}->[$end_index]->[$length] += 1;
		$hits_sum_decimal_len{$chr}->[$end_index]->[$length] += $hit;
		
	}else{
		die $_, $start_index, $end_index;
	}
}
close(IN);

foreach my $chr (sort keys %chr_len){
	my $last_index = int ( ($chr_len{$chr} - 1) / $bin_size);
	for my $i(0.. ($last_index - 1)){
		my $s   = $i * $bin_size + 1;
		my $e	= $s + $bin_size - 1;
		my $cor = "$chr:$s-$e";
		my ($reads, $hits, $HNA) = ( 0, 0, 0 );
		
		if(defined  $raw_reads_sum_integer{$chr}->[$i]){ $reads = $raw_reads_sum_integer{$chr}->[$i] 	}
		if(defined $hits_sum_decimal{$chr}->[$i] ) { $hits =  $hits_sum_decimal{$chr}->[$i]}
		$HNA = $hits * $factor;
	#	print OUT join("\t", ($cor, $reads, $hits, $HNA)), "\n";
		print OUT join("\t", ($cor, $reads, $hits, $HNA)), "\t";
		
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
		
		die $cor if ($total != $reads);
		
	}
	
	my $s   = $last_index * $bin_size + 1;
	my $e	= $chr_len{$chr} ;
	my $cor = "$chr:$s-$e";
	
	my ($reads, $hits, $HNA) = ( 0, 0, 0 );
#print OUT join("\t", ( "coordinate", "raw_reads_sum_integer", "hits_sum_decimal", "HNA"  ) ), "\n";
	
	if(defined  $raw_reads_sum_integer{$chr}->[$last_index]){ $reads = $raw_reads_sum_integer{$chr}->[$last_index] 	}
	if(defined $hits_sum_decimal{$chr}->[$last_index] ) { $hits =  $hits_sum_decimal{$chr}->[$last_index]}
	$HNA = $hits * $factor;
	#print OUT join("\t", ($cor, $reads, $hits, $HNA)), "\n";
	print OUT join("\t", ($cor, $reads, $hits, $HNA)), "\t";
	
	my $i = $last_index;
	
	for my $len ( $bp_interval[0] .. $bp_interval[ $#bp_interval - 1 ] ){
		my ( $reads_sub, $hits_sub, $HNA_sub ) = (0, 0);
		if(defined  $raw_reads_sum_integer_len{$chr}->[$i]->[$len]  ){ $reads_sub = $raw_reads_sum_integer_len{$chr}->[$i]->[$len] 	}
		if(defined  $hits_sum_decimal_len{$chr}->[$i]->[$len]       ){ $hits_sub  = $hits_sum_decimal_len{$chr}->[$i]->[$len]       }
		$HNA_sub = $hits_sub * $factor;
		print OUT join("\t", ( $reads_sub, $HNA_sub )), "\t";
	}
	my $len = $bp_interval[ $#bp_interval];
	my ( $reads_sub, $hits_sub, $HNA_sub ) = (0, 0);
	if(defined  $raw_reads_sum_integer_len{$chr}->[$i]->[$len]  ){ $reads_sub = $raw_reads_sum_integer_len{$chr}->[$i]->[$len] 	}
	if(defined  $hits_sum_decimal_len{$chr}->[$i]->[$len]       ){ $hits_sub  = $hits_sum_decimal_len{$chr}->[$i]->[$len]       }
	$HNA_sub = $hits_sub * $factor;
	print OUT join("\t", ( $reads_sub, $HNA_sub )), "\n";
}


close(OUT);

exit;

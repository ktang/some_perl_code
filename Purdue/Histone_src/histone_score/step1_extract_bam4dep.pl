#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;
my $print_debug = 0;

if ( !$debug ) {
	$print_debug = 0;
}


#my $half_interval_length = 2000;
my $half_flanking_bin_num = 100; #40;
my $bin_size = 50;


if($debug){
	print STDERR "debug = 1\n\n";
}
#my $usage = "$0 \n <input> <output>\n\n";
#die $usage unless(@ARGV == 2);

my $usage = "$0 \n <input_list_DMR> <bam_dir> <outdir> <output_pre>\n\n";
die $usage unless(@ARGV == 4);


my $input	= shift or die;
my $bam_dir	= shift or die;
my $outdir	= shift or die;
my $outpre	= shift or die;

#my $output	= shift or die;
#my $output = File::Spec->catfile($outdir, $outpre . "_half" . $half_interval_length . "bp" . "_bin" . $bin_size ."bp.txt");
my $output = File::Spec->catfile($outdir, $outpre . "_halfBinN" . $half_flanking_bin_num  . "_bin" . $bin_size ."bp.txt");

die unless (-e $input);
die if( -e $output);
die unless (-d $bam_dir);
die if(-e $output);

###########
# get bam file ans label
#############
my @bam_files_simple;
my @bam_files_full;

opendir (DIR, $bam_dir);
@bam_files_simple = grep /\.bam$/, readdir DIR;
closedir DIR;

@bam_files_full = map { File::Spec->catfile($bam_dir, $_ )} @bam_files_simple;

my $file_num = $#bam_files_simple + 1;

foreach my $tmp ( @bam_files_full ){
	die $tmp unless (-e $tmp);
}
open(OUT, ">$output") or die "cannot open $output: $!";

print OUT "#";
print OUT join("\t",@bam_files_simple ), "\n";

if ($debug) {
	print STDERR join("\n", @bam_files_simple ),"\n\n";
	print STDERR join("\n", @bam_files_full ),"\n\n";
	exit;#code
}


###########
# 1, read_list_and_get_interval
##############
my @DMRs_list;# push @DMRs_list ( @{$ref_a}), [$interval_label, $interval_extend ]
read_list ($input, \@DMRs_list,  $bin_size, $half_flanking_bin_num);

#Chr1	3	2

for my $i (0..$#DMRs_list){
	my ($interval_label, $interval_extend ) = @{$DMRs_list[$i]};
	my ($extend_start, $extend_end);
	if ($interval_extend =~ /:(\d+)-(\d+)/) {
		 ($extend_start, $extend_end) = ($1, $2);
	}
	
	#if ($print_debug) { print STDERR join("\t", ()), "\n\n"; }
	if ($print_debug) {
		print STDERR "label1: ", join("\t", ($interval_label, $interval_extend, $extend_start, $extend_end)), "\n" ;#code
	}
	
	my $cmd = "samtools depth -r $interval_extend @bam_files_full";
	my @deps = `$cmd`;
	
	
	my %dep_num = ();
	foreach my $line (@deps){
		#if ($print_debug) { print STDERR  "line: $line;;"  }

		chomp $line;
		my @a = split "\t", $line;
		my $chr = simple_chr($a[0]);
		my $pos = $a[1];
		
		if ($print_debug) { print STDERR join("\t", ($chr, $pos,$extend_start, @a)), "\n";   }
		
		my $index = int ( ($pos - $extend_start) / $bin_size) + 1;
		if ($print_debug) { print STDERR  "index: $index" , "\n"; }
		for my $j (2..$#a){
			$dep_num{$index}->[$j-2] += $a[$j];
			if ($print_debug) {print STDERR $dep_num{$index}->[$j-2] , "\t" };
		}
		if ($print_debug) { print STDERR  ";\n\n"; }

	}
	
	print OUT ">", $interval_label, "\n";
	for my $j (1..($half_flanking_bin_num * 2 +1)){
		my $curr_start = $extend_start + $bin_size * ($j - 1);
		my $curr_end = $curr_start + $bin_size - 1;
		if (! defined $dep_num{$j} ) {
			for my $tmp (0..$#bam_files_simple){
				push @{$dep_num{$j}}, 0;
			}
		}
		
		print OUT join("\t", ($j, $curr_start . "-" .  $curr_end, @{$dep_num{$j}} )), "\n";
	}
	
}

close(OUT);

exit;

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
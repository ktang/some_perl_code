#!/usr/bin/perl -w
## give several mapped bam files, count number for each siRNA seq

use strict;
use File::Spec;

my $debug = 0;
if ($debug) {
	print STDERR "debug = 1 \n\n";
}


# delete sum number < $min_sum
my $min_sum = 4;

my $usage = "$0 \n <outdir> <output_pre> <bam file1> <label1> [<bam2> <label2>.....] \n\n ";
die $usage unless(@ARGV >= 4);

my $outdir = shift or die;
my $outpre = shift or die;
die "outdir not exist" unless (-d $outdir);

my $output = File::Spec->catfile($outdir, $outpre . "_minS_" . $min_sum . ".txt");
die if (-e $output);

my @bams;
my @labels;

my $args_last_index = @ARGV;
my $file_num = int( ( $args_last_index + 1 )/ 2 );

if ($debug) {
	print STDERR join("\n", ($args_last_index, $file_num)), "\n"
}

for my $i (0..($file_num - 1)){
	my $tmp_bam = shift or die;
	my $tmp_label = shift or die;
	die unless (-e $tmp_bam);
	
	push @bams,  $tmp_bam;
	push @labels, $tmp_label; 
}

print STDERR  "bams", "\n";
print STDERR  join("\n", @bams), "\n\n";

print STDERR  "labels", "\n";
print STDERR  join("\n", @labels), "\n";

if ($debug) {
	exit;#code
}


my %siRNAs;
my %siRNA_pos;

open(OUT, ">>$output") or die;
print OUT join("\t", ("seq", "pos", @labels)), "\n";

#FCC159LACXX:6:2203:7745:16899#0/1       16      1       10      255     24M     *       0       0       CTAAACCCTAAACCCTAAACCTCT        IIIIIIIIIIIIIIIIIIIIIIII        XA:i:0  MD:Z:24 NM:i:0  NH:i:1
#col 9
foreach my $i (0..$#bams){
	my $input = $bams[$i];
	my $label = $labels[$i];
	open( IN, "samtools view $input |") or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $seq = $a[9];
		$siRNAs{$seq}->{$label} ++;
		my $pos = join(":", @a[2..3]);
		if (defined  $siRNA_pos{$seq}) {
			die $_ if ($siRNA_pos{$seq} ne $pos);
		}else{
			$siRNA_pos{$seq} = $pos;
		}
		
	}
}

foreach my $k(keys %siRNAs){
	my @tmp_num = ();
	foreach my $l (@labels){
		unless (defined $siRNAs{$k}->{$l}){
			$siRNAs{$k}->{$l} = 0;
		}
		push @tmp_num, $siRNAs{$k}->{$l};
	}
	my $total = get_sum(\@tmp_num);
	if ($total < $min_sum) {
		delete $siRNAs{$k};
		delete $siRNA_pos{$k};
	}
}

foreach my $k(sort keys %siRNAs){
	my @tmp_num = ();
	foreach my $l (@labels){
		push @tmp_num, $siRNAs{$k}->{$l};
	}
	print OUT join("\t", ($k, $siRNA_pos{$k}, @tmp_num)), "\n";

}


exit;

#my $total = get_sum(\@tmp_num);

sub get_sum{
	my ($ref) = @_;
	my $s = 0;
	my $n = scalar(@{$ref});
	for my $i(0..($n-1)){
		$s+=$ref->[$i];
	}
	return $s;
}
#!/usr/bin/perl -w
# two bed files as input
# foreach region in first file, check whether there is a region in overlap with it
# output

# v1.2 output the list

#v1.3 only output the one list, first file overlap with the second
use strict;

my $debug = 0;
print STDERR "debug: $debug\n";

my $usage = "$0 <indir> <bed_file1> <bed_file2> <gap> STDOUT";
#my $usage = "$0 <indir> <bed_file1> <bed_file2> <gap> <outdir> <out_pre>";

die $usage unless (@ARGV == 4);
#die $usage unless (@ARGV == 6);

my ($indir, $input1, $input2, $allowed_gap) = @ARGV[0..3];

#my ($outdir, $f1_over, $f2_over, $f1_left, $f2_left) = @ARGV[4..8];
#my ($outdir, $out_pre) = @ARGV[4..5];

#$ = $out_pre."";

#die "wrong directory" unless (-d $indir and -d $outdir);
die "wrong directory" unless (-d $indir );
#die "wrong input" unless (-e "$indir/$input1" and -e "$indir/$input2" );

#die "wrong output" unless ( !(-e "$outdir/$f1_over") and !(-e "$outdir/$f2_over") and !(-e "$outdir/$f1_left") and !(-e "$outdir/$f2_left") );

open (IN1, "$indir/$input1") or die "cannot open $input1: $!";
open (IN2, "$indir/$input2") or die "cannot open $input2: $!";
#open (OUT, ">$output") or die "cannot open $output: $!";

#open ( O1, ">$outdir/$f1_over") or die "cannot open $f1_over :$!";
#open ( O2, ">$outdir/$f2_over") or die "cannot open $f2_over :$!";
#open ( L1, ">$outdir/$f1_left") or die "cannot open $f1_left :$!";
#open ( L2, ">$outdir/$f2_left") or die "cannot open $f2_left :$!";

my %records;
my %beds;
my $num_f2 = 0;

#my ($head1, $head2);

my $h2 = <IN2>;

while (<IN2>){
#	next if (/start/i);
	$num_f2++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1] + 0;
	my $end = $a[2] + 0;
	
	$records{$chr}->{$start} = $_;
	
	for my $i (($start - $allowed_gap)..($end + $allowed_gap)){
		$beds{$chr}->[$i] = 1;#[$chr, $start];
	} 
}

my $num_f1 = 0;
my $num_overlap = 0;
my $h1 = <IN1>;
print $h1;
while (<IN1>){
#	next if (/start/i);
	$num_f1 ++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $flag = 0;
	#my ($former_chr, $former_start);
	for my $i ($start..$end){
		if(defined $beds{$chr}->[$i]){
			$flag = 1;
			#($former_chr, $former_start) = @{$beds{$chr}->[$i]};
			last;
		}
	}
	if ($flag == 1){
		$num_overlap ++;
		print $_, "\n";
	}
}

my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);

print STDERR "$input2: $num_f2\n";
print STDERR "$input1: $num_f1\noverlap (gap:$allowed_gap) with ", $input2, ": $num_overlap ($per", "%", ")\n";

exit;
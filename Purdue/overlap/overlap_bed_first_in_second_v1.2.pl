#!/usr/bin/perl -w
# two bed files as input
# foreach region in first file, check whether there is a region in overlap with it
# output

# v1.2 output the three lists
# f1 overlapped with f2
# f1 only
# f2 only
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug: $debug\n\n\n";
}

my $usage = "$0 \n <bed_file1> <label1> <bed_file2> <label2> <outdir> \n\n";
#my $usage = "$0 <indir> <bed_file1> <bed_file2> <gap> <outdir> <out_pre>";

#die $usage unless (@ARGV == 9);
die $usage unless (@ARGV == 5);

my $input1 = shift or die;
my $label1 = shift or die;
my $input2 = shift or die;
my $label2 = shift or die;

my $outdir = shift or die;

die $input1 unless (-e $input1);
die $input2 unless (-e $input2);

die "wrong outdir" unless (-d $outdir);

my $output1_only = File::Spec->catfile($outdir, $label1 . "_only_not_overlap_" . $label2 . ".txt");
my $output2_only = File::Spec->catfile($outdir, $label2 . "_only_not_overlap_" . $label1 . ".txt");
my $overlap_file = File::Spec->catfile($outdir, $label1 . "_overlap_" . $label2 . ".txt");

die $output1_only if (-e $output1_only);
die $output2_only if (-e $output2_only);
die $overlap_file if (-e $overlap_file );

if($debug){
	print STDERR join("\n", ( $output1_only, $output2_only, $overlap_file )) , "\n\n";
	exit;
}


open(O1, ">>$output1_only") or die;
die if (-e $overlap_file);
open(OUT, ">>$overlap_file") or die;


open (IN2, $input2) or die $input2;
my %records;
my %coors;
my $num_f2 = 0;

#my ($head1, $head2);
my $h2 = <IN2>;

while (<IN2>){
	
	$num_f2++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1] + 0;
	my $end = $a[2] + 0;
	
	#$records{$chr}->{$start} = $_;
	$coors{$chr}->{$start} = $end;
 		 
	for my $i ($start..$end){
		#$beds{$chr}->[$i] = [$chr, $start];
		$records{$chr}->[$i] = $start;
	}  
}
close(IN2);

open (IN1, $input1 ) or die $input1;
my $num_f1 = 0;
my $num_overlap = 0;

my $head = <IN1>;
print O1 $head;
chomp $head;
my $sym = "overlap_". $label2;

print OUT $head, "\t", $sym, "\n";

while (<IN1>){
	
	$num_f1 ++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $flag = 0;
	#my ($former_chr, $former_start);
	my ($overlap_start, $overlap_end);
	for my $i ($start..$end){
		if(defined $records{$chr}->[$i]){
			$flag = 1;
			$overlap_start = $records{$chr}->[$i];
			$overlap_end  = $coors{$chr}->{$overlap_start};
		#	($former_chr, $former_start) = @{$beds{$chr}->[$i]};
			last;
		}
	}
	
	if ($flag == 1){
		$num_overlap++;
		print OUT join("\t",(@a, join("_", ($chr,$overlap_start, $overlap_end ))) ) , "\n";
	}
	else{
		print O1 $_, "\n";
	}
}
close(IN1);
close(OUT);
close(O1);


my $per = sprintf ("%.1f",100 * $num_overlap / $num_f1);
print STDERR "$label2: $num_f2\n";
print STDERR "$label1: $num_f1\noverlap  with ", $label2, ": $num_overlap ($per", "%", ")\n\n\n";

die  $output2_only if(-e $output2_only);
open(L2, ">> $output2_only") or die;

open (IN1, $input1) or die $input1;
%records = ();
%coors = ();
$num_f1 = 0;

my $h1 = <IN1>;

while (<IN1>){
	
	$num_f1++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1] + 0;
	my $end = $a[2] + 0;
	
	#$records{$chr}->{$start} = $_;
	$coors{$chr}->{$start} = $end;
 		 
	for my $i ($start..$end){
		#$beds{$chr}->[$i] = [$chr, $start];
		$records{$chr}->[$i] = $start;
	}  
}
close(IN1);

open (IN2, $input2 ) or die $input2;
$num_f2 = 0;
$num_overlap = 0;

$head = <IN2>;
print L2 $head;

while (<IN2>){
	
	$num_f2 ++;
	chomp;
	my @a = split /\t/ ;
	my $chr = lc $a[0];
	my $start = $a[1];
	my $end = $a[2];
	my $flag = 0;
	#my ($former_chr, $former_start);
	my ($overlap_start, $overlap_end);
	for my $i ($start..$end){
		if(defined $records{$chr}->[$i]){
			$flag = 1;
			$overlap_start = $records{$chr}->[$i];
			$overlap_end  = $coors{$chr}->{$overlap_start};
		#	($former_chr, $former_start) = @{$beds{$chr}->[$i]};
			last;
		}
	}
	
	if ($flag == 1){
		#print OUT join("\t",(@a, join("_", ($chr,$overlap_start, $overlap_end ))) ) , "\n";
		$num_overlap++;
	}
	else{
		print L2 $_, "\n";
	}
}
close(IN2);
close(L2);

$per = sprintf ("%.1f",100 * $num_overlap / $num_f2);
print STDERR "$label1: $num_f1\n";
print STDERR "$label2: $num_f2\noverlap  with ", $label1, ": $num_overlap ($per", "%", ")\n\n\n";

exit;
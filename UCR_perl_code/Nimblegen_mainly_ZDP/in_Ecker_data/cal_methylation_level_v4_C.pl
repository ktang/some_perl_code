#!/usr/bin/perl -w

#the algorithm is first to calculate the number in the short window( move_size),
#then calculate the number in long window($wsize).

#notice: $wsize/$move_size should be an int.

use strict;

my $usage = "$0 <my_CG_count_chr_data> <Ecker_data> <output_file> <window_size> <move_seze>";

die $usage unless(@ARGV == 5);

open (CG, "<$ARGV[0]")
	or die "cannot open $ARGV[0]:$!";
	
open (ECKER, "<$ARGV[1]")
	or die "cannot open $ARGV[1]:$!";

open (OUT,  ">$ARGV[2]")
  or die "Cannot open $ARGV[2]: $!";

my $wsize = $ARGV[3];
my $move_size = $ARGV[4];
my @pts;
my ($line,$chr,$position,$idx,$index);
my (%nums_move_size, %nums_wsize);



while ( $line = <ECKER>){
	$i_debug++;
	chomp $line;
	@pts = split "\t",$line;
	$position = $pts[3];
	$chr = $pts[0];
	if (($pts[6] > 0)){
		if ( ($position % $move_size) == 0){
			$idx = $position/$move_size - 1;
		}
		elsif ( $position <= 0){
			$idx = 0;
		}
		else{
			$idx = int ($position/$move_size );
		}
		$index = $idx * $move_size + 1; 
		$nums_move_size{$chr}->{$index}++;
		
	}
}

my ($start,$end);

while ($line = <CG>){
	chomp $line;	
	@pts = split "\t",$line;
	$chr = $pts[0];
	$start = $pts[3];
	$end = $pts[4];
	for (my $i = $start;$i <= $end; $i += $move_size){
		if ( (exists $nums_move_size{$chr}) and (exists $nums_move_size{$chr}->{$i})){
			$nums_wsize{$chr}->{$start} += $nums_move_size{$chr}->{$i};
		}
	}
}
close (CG);

my $num_CG;

open (CG, "<$ARGV[0]")
	or die "cannot open $ARGV[0]:$!";
while ($line = <CG>){
	chomp $line;	
	@pts = split "\t",$line;
	$chr = $pts[0];
	$start = $pts[3];
	$end = $pts[4];
	if ( (exists $nums_wsize{$chr}) and (exists $nums_wsize{$chr}->{$start})){
		$num_CG = $nums_wsize{$chr}->{$start};
	}
	else{
		$num_CG = 0;
	}
	$pts[6] = $num_CG;
	if ( $pts[5] == 0){
		if ($num_CG == 0){
			$pts[7] = 0;
		}
		else{
			print "ERROR:$line\nCG_num = $num_CG\n";
		}
	}
	else{
		$pts[7] = $num_CG/$pts[5];
		if ( $pts[7] > 1){
			print "ERROR: percentage > 1\n $line\n$num_CG\n";
			$pts[7] = 1;
		}
	}
	my $out_line = join("\t",@pts);
	print OUT "$out_line\n";
}
close (CG);

close (OUT);
print STDERR "\a";
exit;
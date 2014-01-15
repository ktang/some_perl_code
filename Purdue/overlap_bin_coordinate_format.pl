#!/usr/bin/perl -w
# 1:1-100

use strict;
use File::Spec;

my $debug = 0;
if ($debug){
	print STDERR "debug: $debug\n";
}

#print STDERR "input must have a head\n\n";
my $usage = "$0 \n<input1> <input2> \n\n";
die $usage unless (@ARGV == 2);

my $input1 = shift or die;
my $input2 = shift or die;

die unless (-e $input1);
die unless (-e $input2);

my %records;

my $num_f1 = read_file1 ( $input1, \%records );
my ($num_f2, $num_overlap) = read_file2 ($input2, \%records);

my $per1 = sprintf ("%.1f",100 * $num_overlap / $num_f1);
my $per2 = sprintf ("%.1f",100 * $num_overlap / $num_f2);

# ($volume,$directories,$file) = File::Spec->splitpath( $path );

my ($volume1,$directories1,$file1) = File::Spec->splitpath( $input1 );
my ($volume2,$directories2,$file2) = File::Spec->splitpath( $input2 );

print STDERR  "$file1: $num_f1\t$per1%\n";
print STDERR "$file2: $num_f2\t$per2%\n";
print STDERR  "overlap: ", $num_overlap , "\n\n";
exit;

sub read_file1{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	my $i = 0;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my $coor = $a[0];
		die $_ if( defined $ref->{$coor} );
		 $ref->{$coor} = 1;
	}
	
	close IN;
	return ($i);
}

sub read_file2{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	my $total = 0;
	my $i = 0;

	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		my $coor = $a[0];
		if( defined $ref->{$coor} ){
			$total ++;
		}
		$ref->{$coor} = 1;
	}
	close IN;
	return ($i, $total);

}
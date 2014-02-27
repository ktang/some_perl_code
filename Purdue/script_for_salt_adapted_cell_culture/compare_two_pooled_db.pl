#!/usr/bin/perl -w


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use strict;

my $usage = "\n$0 <input_bigger_DB> <input_small_DB> STDOUT \n\n";
die $usage unless (@ARGV == 2);

my $input_big = shift or die;
my $input_small = shift or die;

#die unless (-e $input);
die unless (-e $input_big);
die unless (-e $input_small);

my %pos_h;
read_big_db($input_big, \%pos_h);

print join("\t", ("chr", "pos")), "\n";


open(IN, $input_small) or die;
while (<IN>) {
	chomp;
	my @a = split "\t";
	my ($chr, $pos) = @a[0..1];
	$chr = Kai_Module::simple_chr($chr);
	unless ( defined $pos_h{$chr}->{$pos} ) {
		print join("\t", ($chr, $pos)), "\n";
	}
	
}

close IN;

exit;

#read_big_db($input_big, \%pos_h);
sub read_big_db{
	my ($file, $ref) = @_;
	die unless ( -e $file );
	
	open(IN, $file) or die;
	my $h = <IN>;
	while (<IN>) {
		chomp;
		my @a = split;
		my ($chr, $pos) = @a[0..1];
		$chr = Kai_Module::simple_chr($chr);
		$ref->{$chr}->{$pos} = 1;
	}
	close IN;
}
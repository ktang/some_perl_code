#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <input_col> <have_head (y or n)> <label>\n\n";
die $usage unless(@ARGV == 4);

#print STDERR "all two files should have head and col start as 1-based\n\n";

my $input_db = shift or die;
my $col_num_input = shift or die;

my $head = shift or die;
die $usage unless ( $head eq "y" or $head eq "n" ); 

my $label = shift or die;

my $index = $col_num_input - 1;


open(IN, $input_db) or die "cannot open $input_db: $!";
if ($head eq "y") {
	my $tmp = <IN>;#code
}

my %n;

while(<IN>){
	chomp;
	my @a = split "\t";
	$n{$a[$index]}++;
}

close(IN);

print join("\t", ("feature", $label, "%". $label)), "\n";

my $t = 0;
for my $k (sort keys %n){
	$t+= $n{$k};
	#print join("\t", ($k, $n{$k})), "\n";
}

for my $k (sort keys %n){
#	$t+= $n{$k};
	my $p = sprintf("%.2f", 100*$n{$k} / $t);
	print join("\t", ($k, $n{$k}, $p)), "\n";
}

print join("\t", ("total", $t, 100)), "\n";

exit;
#record_ID ($ID_file, $ID_file_col_num, \%IDs );
sub record_ID{
	my ($file, $col_num, $ref) = @_;
	die unless (-e $file);
	my $index =  $col_num - 1;
	open(IN, $file);
	while (<IN>){
		chomp;
		my @a = split "\t";
		$ref->{$a[$index]} = 1;
	}
	close IN;
}


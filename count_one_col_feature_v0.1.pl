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
my $usage = "$0 \n <have_head (y or n)> <input> <label> <input_col> [ <input><label> <input_col> ... ]\n\n";
die $usage unless(@ARGV >= 4);

#print STDERR "all two files should have head and col start as 1-based\n\n";


my $head = shift or die;
die $usage unless ( $head eq "y" or $head eq "n" );

my (@inputs, @indexs, @labels) ;

while (@ARGV != 0) {
	my $input = shift or die;
	my $label = shift or die;
	my $col_num_input = shift or die;
	my $index = $col_num_input - 1;

	
	push  @inputs, $input;
	push  @indexs, $index;
	push  @labels, $label;

}

die unless ( @inputs == @indexs and @inputs == @labels);

if ($debug) {
	for my $i(0..$#inputs){
		print STDERR join("\t", ($inputs[$i], $labels[$i],$indexs[$i])), "\n";
	}
	exit;
}



my %n;

for my $i(0..$#inputs){
	my $input_db = $inputs[$i];
	my $label = $labels[$i];
	my $index = $indexs[$i];
	open(IN, $input_db) or die "cannot open $input_db: $!";
	if ($head eq "y") {
		my $tmp = <IN>;#code
	}


	while(<IN>){
		chomp;
		my @a = split "\t";
		$n{$a[$index]}->{$label}++;
	}

	close(IN);
}


my @per_labels = map {"%" . $_} @labels;

print join("\t", ("feature", @labels, @per_labels)), "\n";

my %t ;
for my $k (sort keys %n){
	for my $l ( @labels){
		unless(defined $n{$k}->{$l}){
			$n{$k}->{$l} = 0;
		}
		$t{$l}+=$n{$k}->{$l};
	}
}

for my $k (sort keys %n){
	my @n_tmp = ();
	my @per_tmp = ();
	for my $l ( @labels){
		my $p = sprintf("%.2f", 100 * $n{$k}->{$l} / $t{$l} );
		push @per_tmp, $p;
		push @n_tmp, $n{$k}->{$l};
	}
	
	print join("\t", ($k, @n_tmp, @per_tmp)), "\n";
}

#print join("\t", ("total", $t, 100)), "\n";

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


#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_db> <input_col_as_ID> <ID_file> <col_as_ID> <output>\n\n";
die $usage unless(@ARGV == 5);

print STDERR "all two files should have head and col start as 1-based\n\n";

my $input_db = shift or die;
my $col_num_input = shift or die;
my $index_db = $col_num_input - 1;

my $ID_file = shift or die;
my $ID_file_col_num = shift or die;
my $output = shift or die;

die unless (-e $input_db);
die if( -e $output and $output ne "/dev/null");

my %IDs;

record_ID ($ID_file, $ID_file_col_num, \%IDs );



open(IN, $input_db) or die "cannot open $input_db: $!";
die if(-e $output and $output ne "/dev/null");
open(OUT, ">$output") or die "cannot open $output: $!";
while(<IN>){
    next if /^#/;
    last;
}
#my $h = <IN>;
#print OUT $h;
print OUT $_;
while(<IN>){
	chomp;
	my @a = split "\t";
	if(defined $IDs{$a[$index_db]}){
		print OUT $_, "\n";
		delete $IDs{$a[$index_db]};
	}
}

close(IN);
close(OUT);

for my $k (sort keys %IDs){
	print STDERR $k, "\n";
}

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


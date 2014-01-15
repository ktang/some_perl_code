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
#my $usage = "$0 \n <input_db_for_extract> <ID_file_as_target>  <output>\n\n";
my $usage = "$0 \n<ID_file_as_target>  <input_db_for_extract>  <output>\n\n";
die $usage unless(@ARGV == 3);

#print STDERR "all two files should have head and col start as 1-based\n\n";

my $ID_file = shift or die;
my $input_db = shift or die;
my $output = shift or die;

die unless (-e $input_db);
die if( -e $output and $output ne "/dev/null");

my %IDs;

#record_ID ($ID_file, $ID_file_col_num, \%IDs );
record_bed($ID_file,  \%IDs );


open(IN, $input_db) or die "cannot open $input_db: $!";
die if(-e $output and $output ne "/dev/null");
open(OUT, ">$output") or die "cannot open $output: $!";

#while(<IN>){
 #   next if /^#/;
  #  last;
#}
my $h = <IN>;
print OUT $h;
#print OUT $_;

while(<IN>){
	if(/^#/){
		print OUT $_;
		next;
	}
	chomp;
	my @a = split "\t";
	if ($a[0] !~ /chr/) {
		$a[0] = "chr" . $a[0];
	}
	
	my $index = join("_",@a[0..2]);
	if(defined $IDs{$index}){
		print OUT $_, "\n";
		delete $IDs{$index};
	}
}

close(IN);
close(OUT);

for my $k (sort keys %IDs){
	print STDERR $k, "\n";
}

exit;

# record_bed ($ID_file,  \%IDs );
sub record_bed{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(IN, $file);
	while (<IN>){
		chomp;
		my @a = split "\t";
		my $index = join("_",@a[0..2]);
		$ref->{$index} = 1;
	}
	close IN;
}

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


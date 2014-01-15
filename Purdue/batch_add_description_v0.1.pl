#!/usr/bin/perl -w

#v0.1
# input can be AT1G00001.X
#and can have no head

use strict;
use File::Spec;

#my $indir = ".";

my $usage = "$0 \n <indir> <outdir>  <gene_col_number> <head:y or n>\n\n";
die $usage unless (@ARGV == 4);

my $indir = shift or die ;
my $outdir = shift or die;

my $gene_col_number = shift or die;

my $gene_index = $gene_col_number - 1;

my $head_flag = shift or die;

die $usage unless ($head_flag eq "y" or $head_flag eq "n");

opendir(DIR, $indir) or die "cannot open $indir";
my @inputs = grep /.txt$/, readdir DIR;
closedir DIR;

print STDERR join ("\n", @inputs) , "\n\n";

my $func = "/Users/tang58/DataBase/TAIR10/TAIR10_functional_descriptions.txt";

#my $gene_index = 5;

my %func;
open(FF, $func) or die "Can't open $func: $!";
while(<FF>){
	next if(/^Model_name/);
	chomp;
	my @temp = split /\t/;
	if($temp[0] =~ /(AT.G\d+)/){
		my $gene = $1;
		if(!defined $func{$gene}){
			if(defined $temp[2] && $temp[2] =~ /\w+/){
			   $func{$gene} = $temp[2];
			}elsif(defined $temp[3] && $temp[3] =~ /\w+/){
				$func{$gene} = $temp[3];
			}elsif(defined $temp[4] && $temp[4] =~ /\w+/){
				$func{$gene} = $temp[4];
			}else{
				$func{$gene} = "NONE";
			}
		}
	}else{
		#die "gene name ", $temp[0], " does not match pattern";
		last;
	}
}




foreach my $input(@inputs){
	if($input =~ /(\S+)\.txt$/){
		my $output = File::Spec->catfile($outdir, $1."_func.txt");
		
		my $full_input =  File::Spec->catfile( $indir, $input);
		
		open (IN, $full_input) or die;
		
		if (-e $output){
			print STDERR "$output exists\n";
			next;
		}
		
		else{
			open (OUT, ">$output");
			
			if($head_flag eq "y" ){
				my $head = <IN>;
				chomp $head;
				print OUT join("\t", ($head, "func") ), "\n";
			}			
			while (<IN>){
				chomp;
				my @a = split "\t";
				my $gene = "NA";
				if ( $a[$gene_index] =~ /^(AT.G\d\d\d\d\d)$/ ){
					$gene = $1;
				}elsif( $a[$gene_index] =~ /^(AT.G\d\d\d\d\d)\.\d/){
					$gene = $1;
				}else{
					die $a[$gene_index];
				}
				
				if (defined $func{$gene}){
					print OUT join("\t", ($_, $func{$gene} )), "\n";
				}
				else{
					print STDERR $input, "\n";
					print STDERR $_, "\n";
					print STDERR "$a[$gene_index] has no description\n";
				}
			}
			close(OUT);
		}#else
		
		close(IN);
	}
}
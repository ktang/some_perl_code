#!/usr/bin/perl -w

use strict;
use File::Spec;
my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> STDOUT\n\n";
die $usage unless(@ARGV == 1);


my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";

die unless (-d $indir);
#die unless (-d $outdir);

opendir (DIR, $indir) or die;
my @files = grep /\.txt$/ , readdir DIR;
closedir DIR;

print STDERR "files: \n";
print STDERR join("\n", @files ), "\n\n";

my @labels= ();
foreach my $file (@files){
	my $pre;
	
	if ($file =~ /(\S+)_SelfPer/){
		$pre = $1;
		push @labels, $pre;
	}else{
		die $file;
	}
}

my %nums; #h{mut_type}->file = 
#chr	pos	ref	mut_nrpe1_150	dep_nrpe1_150	per_nrpe1_150	seq_nrpe1_150	qual_nrpe1_150
#chr1	29081	G	G=>DEL	14	14.29	...,,,.,,,,-2ta,-2ta.,	DFFI)@J>##D4H<

my $col_num = 4;

if ($debug) {
	exit;#code
}


foreach my $file (@files){
	my $pre;
	
	if ($file =~ /(\S+)\.txt$/){
		$pre = $1;
	}else{
		die $file;
	}
	
	my $input = File::Spec->catfile($indir, $file);
	
	open(IN, $input) or die "$input: $!";
	my $head = <IN>;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $type = $a[ $col_num  - 1];
		$nums{$type}->{$file} += 1;
	}
	
	close IN;
}

print join("\t", ("mutation_type", @labels)), "\n";
foreach my $type (sort keys %nums){
	print $type, "\t";
	foreach my $file(@files){
		my $n = 0;
		if (defined $nums{$type}->{$file}) {
			$n = $nums{$type}->{$file};#code
		}
		print $n;
		if ($file =~/^6/) {
			print "\n";
		}else{
			print "\t";
		}
		
	}
}
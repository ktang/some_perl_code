#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
print STDERR "debug:$debug\n";

my $usage = "$0 <indir> <outdri> <chr|Chr>";

die $usage unless (@ARGV == 3);

#my ($indir, $flag) = @ARGV[0..1];
my $indir = shift or die;
my $outdir = shift or die;
my $flag = shift or die;

die "wrong indir" unless (-d $indir);
die unless (-d $outdir);

die "Chr|chr" unless ($flag eq "Chr" or $flag eq "chr");

my $C24_file = "";
if($flag eq "Chr"){
	$C24_file = "/Users/tang58/deep_seq_analysis/SNP_pos_in_C24andRos1_Chr.txt";
}
elsif ($flag eq "chr"){
	$C24_file = "/Users/tang58/deep_seq_analysis/SNP_pos_in_C24andRos1.txt";
}
else{
	die "wrong flag";
}

#if ($debug){
#	$C24_file = "/Users/tang58/try/debug/C24_SNP_debug_Chr.txt";
#}

opendir (DIR, $indir) or die "cannot open $indir";
my @files = grep /pileup$/, readdir  DIR;

print STDERR join("\n", @files), "\n";

open (C24, $C24_file) or die "cannot open C24:$!";
my %SNPs;

foreach my $infile (@files){
	if($infile =~ /(\S+)\.pileup$/){
		my $pre = $1;
		
		my $output = File::Spec->catfile($outdir, $pre."_not_C24.txt");
		my $input =  File::Spec->catfile($indir, $infile);
		
		die "wrong input" unless (-e $input );
		die "wrong output" if (-e $output );
		
		if($debug){
			print join("\t", ($input, $output)), "\n\n";
		}
		
	}
	else{
		print STDERR "wrong input $infile\n";
		die;
	}	
}

if($debug){
	die "debut";
}

print STDERR "reading $C24_file...\t";
while(<C24>){
	chomp;
	my @a = split "\t";
	my $chr = $a[0];
	my $pos = $a[1];
	$SNPs{$chr}->[$pos] = 1;
}
print STDERR "Done\n";

foreach my $infile (@files){
	
	print STDERR "handling $infile...\t";
	
	if($infile =~ /(\S+)\.pileup$/){
		my $pre = $1;
		my $output = File::Spec->catfile($outdir, $pre."_not_C24.txt");
		my $input =  File::Spec->catfile($indir, $infile);
		
		die "wrong input" unless (-e $input );
		die "wrong output" if (-e $output );
				
		open (IN, $input) or die "cannot open $input:$!";
		open (OUT, ">$output") or die "cannot open $output:$!";

		while(<IN>){
#			next if /^\#/;
#			chomp;
			my @a = split "\t";
			my $chr = $a[0];
			my $pos = $a[1];
			print OUT $_ unless (defined $SNPs{$chr}->[$pos]);
		}

		close(IN);
		close(OUT);
	}
	else{
		print STDERR "wrong input $infile\n";
	}	
	print STDERR "Done\n";
}

exit;
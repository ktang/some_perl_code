#!/usr/bin/perl -w
#v0.1: use File::Spec; and add outdir

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug:$debug\n";
}
print STDERR "Chr\n";

#my $usage = "$0 <indir> <chr|Chr>";
#die $usage unless (@ARGV == 2);
#my ($indir, $flag) = @ARGV[0..1];

my $usage = "$0 <indir> <outdir>";
die $usage unless (@ARGV == 2);
my $indir = shift or die "indir";
my $outdir = shift or die "outdir";


die "wrong indir" unless (-d $indir);

#die "Chr|chr" unless ($flag eq "Chr" or $flag eq "chr");


my $exon_gff = "/Users/tang58/DataBase/TAIR10/GFF/TAIR10_GFF3_exon.gff";

opendir (DIR, $indir) or die "cannot open $indir";
my @files = grep /(\S+)_not_C24\S+\.txt$/, readdir  DIR;

print STDERR join("\n", @files), "\n";



foreach my $file (@files){
	if($file =~ /(\S+)_not_C24\S+\.txt$/){
		my $pre = $1;
		my $output = File::Spec->catfile($outdir, $pre."_in_exon.txt");
		my $input = File::Spec->catfile($indir, $file);
		
		die "wrong input" unless (-e $input );
		die "wrong output" if (-e $output);
		
		print STDERR join("\t", ($input, $output)), "\n";
	}
	else{
		print STDERR "wrong input $file\n";
	}	
}

if($debug) {
	die "debug\n\n";
}

open (GFF, $exon_gff) or die "cannot open $exon_gff:$!";
my %SNPs;

print STDERR "reading $exon_gff...\t";
while(<GFF>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[0];
#	my $pos = $a[1];
	my $start = $a[3];
	my $end = $a[4];
	for my $pos($start..$end){
		$SNPs{$chr}->[$pos] = 1;
	}
}
print STDERR "Done\n";

foreach my $file (@files){
	
	print STDERR "handling $file...\t";
	
	if($file =~ /(\S+)_not_C24\S+\.txt$/){
		my $pre = $1;
		my $output = File::Spec->catfile($outdir, $pre."_in_exon.txt");
		my $input = File::Spec->catfile($indir, $file);
		
		die "wrong input" unless (-e $input );
		die "wrong output" if (-e $output);
		
		open (IN, $input) or die "cannot open $file:$!";
		open (OUT, ">$output") or die "cannot open $output:$!";

		while(<IN>){
#			next if /^\#/;
#			chomp;
			my @a = split "\t";
			my $chr = lc $a[0];
			my $pos = $a[1];
			print OUT $_ if (defined $SNPs{$chr}->[$pos]);
		}

		close(IN);
		close(OUT);
	}
	else{
		print STDERR "wrong input $file\n";
	}	
	print STDERR "Done\n";
}

exit;
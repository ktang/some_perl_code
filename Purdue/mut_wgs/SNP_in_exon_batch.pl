#!/usr/bin/perl -w

use strict;

my $debug = 0;
print STDERR "debug:$debug\n";
print STDERR "Chr\n";

#my $usage = "$0 <indir> <chr|Chr>";
#die $usage unless (@ARGV == 2);
#my ($indir, $flag) = @ARGV[0..1];

my $usage = "$0 <indir> ";
die $usage unless (@ARGV == 1);
my $indir = $ARGV[0];


die "wrong indir" unless (-d $indir);

#die "Chr|chr" unless ($flag eq "Chr" or $flag eq "chr");


my $exon_gff = "/Users/tang58/DataBase/TAIR10/GFF/TAIR10_GFF3_exon.gff";

opendir (DIR, $indir) or die "cannot open $indir";
my @files = grep /not_C24\.txt$/, readdir  DIR;

print STDERR join("\n", @files), "\n";



foreach my $input (@files){
	if($input =~ /(\S+)_not_C24\.txt$/){
		my $pre = $1;
		my $output = $1."_in_exon.txt";
		die "wrong input" unless (-e "$indir/$input" );
		die "wrong output" unless ( !(-e "$indir/$output"));
	}
	else{
		print STDERR "wrong input $input\n";
	}	
}

open (GFF, $exon_gff) or die "cannot open $exon_gff:$!";
my %SNPs;

print STDERR "reading $exon_gff...\t";
while(<GFF>){
	chomp;
	my @a = split "\t";
	my $chr = $a[0];
#	my $pos = $a[1];
	my $start = $a[3];
	my $end = $a[4];
	for my $pos($start..$end){
		$SNPs{$chr}->[$pos] = 1;
	}
}
print STDERR "Done\n";

foreach my $input (@files){
	
	print STDERR "handling $input...\t";
	
	if($input =~ /(\S+)_not_C24\.txt$/){
		my $pre = $1;
		my $output = $1."_in_exon.txt";
		die "wrong input" unless (-e "$indir/$input" );
		die "wrong output" unless ( !(-e "$indir/$output"));
		
		open (IN, "$indir/$input") or die "cannot open $input:$!";
		open (OUT, ">$indir/$output") or die "cannot open $output:$!";

		while(<IN>){
#			next if /^\#/;
#			chomp;
			my @a = split "\t";
			my $chr = $a[0];
			my $pos = $a[1];
			print OUT $_ if (defined $SNPs{$chr}->[$pos]);
		}

		close(IN);
		close(OUT);
	}
	else{
		print STDERR "wrong input $input\n";
	}	
	print STDERR "Done\n";
}

exit;
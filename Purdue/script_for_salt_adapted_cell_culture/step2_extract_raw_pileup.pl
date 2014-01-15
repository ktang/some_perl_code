#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
my $constant = 100; #100

if ($debug) {
	print STDERR "debug = 1\n\n";#code
}


my $usage = "$0 \n <indir_bam> <input_pos_file> <output>\n\n";
die $usage unless (@ARGV == 3);
my $indir = shift or die;
my $input = shift or die;
my $output = shift or die;

my $fa_file = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa";
die unless (-e $fa_file);

die unless ( -d $indir);
die if(-e $output);

die unless (-e $input);


open(OUT, ">>$output") or die; 
#print OUT join("\t",("chr", "pos") ), "\n";
	       
opendir (DIR, $indir) or die ;
my @files = grep /\.bam$/ , readdir DIR;
closedir DIR;


if ($debug) {
	print STDERR join("\n", @files), "\n\n";
}

my @full_files =();
my @simple_labels = ();
my @heads = ();

foreach my $file(@files){
	my $input = File::Spec->catfile($indir, $file);
	my $bai_file = File::Spec->catfile($indir, $file . ".bai");

	die $bai_file unless ( -e $bai_file);
	die $file unless (-e $input);

	push @full_files, $input ; 
	
	if ( $file =~ /(\S+)_bowtie2_/) {
		my $lab = $1;
		push @simple_labels, $lab;
		push @heads, "d_"     . $lab;
		push @heads, "p_"  . $lab;
		push @heads, "q_"    . $lab;
		
	}else{
		die $file;
	}
}

die "num not equal" unless ( @full_files == @simple_labels);


if ($debug) {
	for my $i (0..$#full_files){
	print STDERR join("\t", ($simple_labels[$i],  $full_files[$i]) ), "\n";
	}
	print STDERR "\n\n";

	print STDERR join("\n",@heads), "\n\n";
}

my $all_bam = join(" ", @full_files);
open(IN, $input) or die "cannot open $input";

print OUT join("\t", ("chr", "pos", "ref", @heads)) , "\n";

if ($debug) {
	exit;#code
}


# samtools mpileup  -f /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa -r 1:973010-973023
#N11_S20_XZ7_ddc_150mM_bowtie2_XH_Dec24_sort_rmdup.bam N01_S27_XZ9_Col0_E_bowtie2_XH_Dec24_sort_rmdup.bam  2>/dev/null 

#my $head = <IN>;

my $max_dep = 1000;

my $i = 0;

while (<IN>) {
	
	chomp;
	my ($chr, $p) = split "\t";
	next if ($chr =~ /^chr$/i);
	$i++;
	if ($i % $constant == 0) {
		print STDERR "$i";
		if ( $i % ( 10 * $constant) == 0 ) {
			print STDERR "\n";
		}else{
			print "\t";
		}
		
	}
	
	my $region = "$chr:$p-$p";
	my $cmd = "samtools mpileup -d $max_dep -f $fa_file -r $region  $all_bam 2> /dev/null"; 
#	my $cmd = "samtools mpileup -d $max_dep -f $fa_file -r $region  $all_bam ";
	my $tmp = `$cmd`;
	print OUT $tmp;
}

close IN;
close OUT;
exit;


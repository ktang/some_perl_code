#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/step4_reduce_boundary_v0.2.pl";
die unless (-e $script );
# <input_list> <format (bed or coordinate)> <outdir> <outpre> <bam_file1> [ <bam2> ...] 

if($debug){
	print STDERR "debug: $debug\n";
}

my $usage = "\n <indir> <outdir> <format (bed or coordinate)> <bam_file1> [ <bam2> ...]  \n\n";
die $usage unless (@ARGV >= 4);

#my ($indir, $outdir, $pattern) = @ARGV;
my $indir = shift or die;
my $outdir= shift or die;
my $formart = shift or die;
#my $XX_bp = shift or die;

die unless ($formart eq "bed" or $formart eq "coordinate");

die "wrong dir" unless (-d $indir and -d $outdir);

my @bams;
for my $tmp(0..$#ARGV){
	my $f = shift or die;
	die unless ( -e $f);
	push @bams, $f;
}


opendir (DIR, $indir) or die "cannot open $outdir: $!";
my @files = grep /\.txt$/, readdir DIR;
closedir(DIR);


if ($debug){
	print STDERR join("\n", @files),"\n\n\n";
	#exit;
}

foreach my $input(@files){
	if( $input =~ /(\S+)\.txt$/){
		my $outpre = $1;
		my $full_input = File::Spec->catfile($indir, $input);
		#my $output = File::Spec->catfile( $outdir,  $1 . "_" . $XX_bp . "bp.bam" );
		
		
		die "wrong input: $full_input" unless (-e $full_input);
		#die "wrong output" if  ( -e $output );
		my $cmd = "perl $script $full_input $formart $outdir $outpre @bams";
		# <input_list> <format (bed or coordinate)> <outdir> <outpre> <bam_file1> [ <bam2> ...]
		print STDERR $cmd, "\n\n";

	
		if($debug){
			#print STDERR $cmd, "\n\n";
			#print STDERR $output, "\n";
		}
		else{
			#print STDERR "handling $input...";
			`$cmd`;
		}
	}
	else{
		die $input;
	}
}

exit;

sub extract_bam{
	my ($input_sub, $output_sub, $XX_bp_sub) = @_;
	die unless ( -e $input_sub);
	die if(-e $output_sub);
	
	open( IN, "samtools view -h $input_sub |") or die;
	open( OUT ,  "| samtools view -b -S -  > $output_sub") or die;
	
	while(<IN>){
		if(/^@/){
			print OUT $_;

		}else{
			#chomp;
			my @a = split "\t";
			print OUT $_ if (length ($a[9]) ==  $XX_bp_sub);
		}
	}
	close IN;
	close OUT

}

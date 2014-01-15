#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug: $debug\n";
}

my $usage = "\n<indir> <outdir> <XX_bp>\n\n";
die $usage unless (@ARGV == 3);

#my ($indir, $outdir, $pattern) = @ARGV;
my $indir = shift or die;
my $outdir= shift or die;
my $XX_bp = shift or die;

die "wrong dir" unless (-d $indir and -d $outdir);

opendir (DIR, $indir) or die "cannot open $outdir: $!";
my @files = grep /\.bam$/, readdir DIR;
closedir(DIR);


if ($debug){
	print STDERR join("\n", @files),"\n\n\n";
	#exit;
}

foreach my $input(@files){
	if( $input =~ /(\S+)\.bam$/){
		
		my $full_input = File::Spec->catfile($indir, $input);
		my $output = File::Spec->catfile( $outdir,  $1 . "_" . $XX_bp . "bp.bam" );
		
		
		die "wrong input: $full_input" unless (-e $full_input);
		die "wrong output" if  ( -e $output );
		
		
		
		#my $cmd = "perl $script $full_input > $output";
		
	
		if($debug){
			#print STDERR $cmd, "\n\n";
			print STDERR $output, "\n";
		}
		else{
			print STDERR "handling $input...";
			extract_bam($full_input, $output, $XX_bp);
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

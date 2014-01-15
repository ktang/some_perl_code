#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/step2.1_cal_HNA_for_non_overlaping_interval_v0.0_all_length.pl";#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);
#my  <input_bed_list> <input_bam> <label> <total_reads>  <output>

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir> <input_bam> <label> <total_reads>\n\n";
die $usage unless(@ARGV == 5);

my $indir = shift or die;
my $outdir = shift or die;

my $inbam = shift or die;
my $label = shift or die;
my $total_reads = shift or die;

die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir) or die "cannot open $indir: $!";
my @files = grep /\.txt$/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	if ($file =~ /(\S+)\.txt$/) {
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		my $output = File::Spec->catfile ($outdir, $pre . "_sRNA.txt");
		my $cmd = "perl $script $input $inbam $label $total_reads $output";
		print STDERR $cmd, "\n\n";
		if (!$debug) {
			`$cmd`;
		}
		
	}
	
}

exit;


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

print STDERR "\n\n please change ada_seq\n\n\n";

my $debug = 0;

#my $ada_seq = "AGATCGGAAGAGCACACGTC";

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/bowtie_pipeline/bowtie1_sRNA_pipeline_Col0_v0.1.pl";
#"/Users/tang58/Kai_BS/PCA_analysis/1_cal_methylation_level_in_bins.pl";
die unless (-e $script);

if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;

die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir) or die "cannot open $indir: $!";
my @files = grep /L001_clean/, readdir DIR;
closedir DIR;

foreach my $file (@files){
	if ( $file =~ /(\S+)_L00/) {
		my $pre = $1 ;
#		my $input = File::Spec->catfile($indir, $file);
#		die unless (-e $input);
		#my $usage = "$0 \n\n <indir> <input_fa_gz_pre> <outdir> <outpre>\n\n";


		my $cmd = "$script $indir $pre  $outdir $pre";
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
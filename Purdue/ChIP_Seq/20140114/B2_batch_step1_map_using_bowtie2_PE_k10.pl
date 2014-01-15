#!/usr/bin/perl -w

use strict;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}
my %phred_label = ( 33=> "--phred33", 64=>"--phred64");

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/ChIP_Seq/20140114/step1_map_using_bowtie2_PE_ChIP_20140114.pl";
die "script" unless (-e $script);

my $usage = "\n$0 \n <indir> <outdir> <phred_33_64> <p> <trim5> <trim3>\n\n";

die $usage unless (@ARGV == 6);
my $indir = shift or die;
die "indir is not a dir" unless (-d $indir);

my $outdir = shift or die;
die "outdir is not a dir" unless(-d $outdir);

my $phred_num = shift or die;
die unless (defined $phred_label{$phred_num});

my $p =      shift or die; #2;
my $trim5 =  shift or die;#3;
my $trim3  = shift or die;#12;


opendir(DIR, $indir) or die "cannot open indir";

#my @files = grep /_L001_R1\.fq\.gz$/, readdir DIR;
my @files = grep /R1\.fq\.gz$/, readdir DIR;

print STDERR join ("\n", @files), "\n\n";

# B1
#      1 35S_SUC2_Input_131216_S03_AAGACG_L001_R1.fq.gz
#      2 35S_SUC2_Input_131216_S03_AAGACG_L001_R2.fq.gz
#      3 35S_SUC2_Input_131216_S03_AAGACG_L002_R1.fq.gz
#      4 35S_SUC2_Input_131216_S03_AAGACG_L002_R2.fq.gz


foreach my $file (@files){
	if ($file =~ /(\S+)_R1\.fq.gz$/){
		my $pre = $1;
		
		my $outpre = "";
		if ( $pre =~ /(\S+)_S\d\d_/) {
			$outpre = $1;
		}else{
			die $pre;
		}
		
		
		my $cmd = "time $script $indir $pre $outdir $outpre $phred_num $p $trim5 $trim3";
		print STDERR $cmd, "\n\n";
		unless($debug){
			`$cmd`;
		}
	}
	else{
		die $file;
	}
}
#my $usage = "\n$0 \n\n<indir> <inpre> <outdir> <outpre> <phred_33_64>\n\n";
exit;
#!/usr/bin/perl -w

# use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/RNA_Seq/cuffdiff_single_sample.pl";
die "script" unless (-e $script);

use strict;
use File::Spec;

my $debug = 0;
my $usage = "$0 \n <WT_bam> <wt_lab> <mut_indir> <all_outdir> <CPU_num>\n\n";
die $usage unless(@ARGV == 5);


my $bam_wt  = shift or die;
my $lab_wt  = shift or die;
my $indir = shift or die "indir";
my $outdir = shift or die;
my $CPU_num = shift or die;

die "indir" unless (-d $indir);
die "outdir" unless(-d $outdir);
die unless (-e $bam_wt);


opendir(DIR, $indir) or die "dir";

my @files = grep /\.bam$/, readdir DIR;

closedir DIR;

print STDERR join("\n", @files), "\n\n";


# my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";

foreach my $file(@files){
	#next if ( $file =~ /$bam_wt/ );
	#next if ($bam_wt =~ /$file/);
	
	if ($bam_wt =~ /$file/ or  $file =~ /$bam_wt/){
		print STDERR $file . "match WT\n\n";
		next;
	}
	
	if($file =~ /(\S+)_thout.+\.bam$/){
		
		my $lab_mut = $1;
		my $sub_dir = File::Spec->catfile($outdir, $lab_mut . "_vs_" . $lab_wt . "_cuffdiff" );
		die if (-d $sub_dir);
		
		my $input = File::Spec->catfile($indir, $file);
		die unless (-e $input);
		
		#my $usage = "$0 \n <mut_bam> <mut_label> <wt_bam> <wt_label> <outdir> <CPU_num>\n\n";

		my $cmd = "perl $script $input $lab_mut $bam_wt $lab_wt  $sub_dir $CPU_num";
		print STDERR $cmd, "\n\n";
		if(!$debug){	
			`$cmd`;
		}else{
			print STDERR "OK\n";
		}
	}
	else{
		die $file;
	}
}

exit;

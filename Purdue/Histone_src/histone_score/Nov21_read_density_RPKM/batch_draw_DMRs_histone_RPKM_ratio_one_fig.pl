#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $R_script = "/Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/Nov21_read_density_RPKM/draw_DMRs_histone_RPKM_ratio_one_fig.r";
die unless (-e $R_script);
my $usage = "$0 \n <indir>    <binSize: 200> <half_bin_num: 25> \n\n";
die $usage unless(@ARGV == 3);

my $indir	= shift or die;
#my $list_pre	= shift or die;
my $binSize	= shift or die;
my $half_bin_num= shift or die;

#die unless (-d $indir);
die unless (-d $indir);


opendir(DIR, $indir) or die "cannot open $indir: $!";
my @DMR_lists = grep /_Input_/ , readdir DIR;
closedir DIR;




if ( $debug ) {
	print STDERR join("\n", @DMR_lists), "\n\n";
	#print STDERR join("\n", @bam_lists), "\n\n";
}

#R --slave --vanilla --args indir list_pre binSize half_bin_num  outpre
foreach my $file(@DMR_lists){
	if ( $file =~ /(\S+)_in_/) {
		my $list_pre = $1;
		my $out_pre = $list_pre . "_ChIP_RPKM_ratio";
		my $cmd = "R --slave --vanilla --args $indir  $list_pre $binSize $half_bin_num  $out_pre < $R_script";
		print STDERR $cmd, "\n\n";
		if ( !$debug) {
			`$cmd`;#code
		}
	}
}
# <input_list_DMR> <bam_dir> <outdir> <output_pre>
print STDERR "\n\n";
exit;

#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 0;

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/Nov21_read_density_RPKM/step1_extract_ChIP_read_num.pl";
die unless (-e $script);
#my $usage = "$0 \n <input_list_DMR> <type: bed or coor> <bam_file>  <bam_label> <read_length> <outdir> <output_pre>\n\n";

if($debug){
	print STDERR "debug = 1\n\n";
}

my %types = ("bed"=>1, "coor" => 1);

my $usage = "$0 \n <DMR_indir>  <type: bed or coor>  <bam_dir> <read_length> <outdir>\n\n";
die $usage unless(@ARGV == 5);

my $indir_DMR = shift or die;
my $type	= shift or die;
die "type: bed or coor" unless (defined $types{$type});
my $bam_dir = shift or die;
my $read_length = shift or die;
my $outdir = shift or die;

#die unless (-d $indir);
die unless (-d $indir_DMR);
die unless (-d $bam_dir);
die unless (-d $outdir);

opendir(DIR, $indir_DMR) or die "cannot open $indir_DMR: $!";
my @DMR_lists = grep /\.txt$/ , readdir DIR;
closedir DIR;

opendir(DIR, $bam_dir) or die "cannot open $bam_dir: $!";
my @bam_lists = grep /\.bam$/ , readdir DIR;
closedir DIR;



if ( $debug ) {
	print STDERR join("\n", @DMR_lists), "\n\n";
	print STDERR join("\n", @bam_lists), "\n\n";
}

#<input_list_DMR> <type: bed or coor> <bam_file>  <bam_label> <read_length> <outdir> <output_pre>
foreach my $file(@DMR_lists){
	if ( $file =~ /(\S+)\.txt$/) {
		my $bed_pre = $1;
		my $bed_list = File::Spec->catfile( $indir_DMR, $file );
		die "$bed_list do not exists" unless (-e $bed_list);
		foreach my $bam_file(@bam_lists){
			my $bam_input = File::Spec->catfile( $bam_dir, $bam_file);
			die "$bam_file do not exists\n\n" unless (-e $bam_input);
			if ( $bam_file =~ /(\S+)\.bam$/) {
				my $bam_pre = $1;
				my $cmd = "perl $script $bed_list $type  $bam_input $bam_pre $read_length $outdir $bed_pre";
				#print STDERR $cmd, "\n\n";
				print STDERR $cmd, "\n\n";
				if ( !$debug) {
					`$cmd`;#code
				}
			}
		}
	}
}
# <input_list_DMR> <bam_dir> <outdir> <output_pre>
print STDERR "\n\n";
exit;

#!/usr/bin/perl -w

# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}


my @chrs = ( "Chr1",  "Chr2",  "Chr3",  "Chr4",  "Chr5", "chloroplast", "mitochondria");

my $usage = "$0 \n <DMR_indir>  <output>\n\n";
die $usage unless(@ARGV == 2);

my $bam_dir = shift or die;
my $output = shift or die;
die "output exists\n" if (-e $output);

opendir(DIR, $bam_dir) or die "cannot open $bam_dir: $!";
my @bam_lists = grep /\.bam$/ , readdir DIR;
closedir DIR;



if ( $debug ) {
	print STDERR join("\n", @bam_lists), "\n\n";
}

#<input_list_DMR> <type: bed or coor> <bam_file>  <bam_label> <read_length> <outdir> <output_pre>

unless($debug){
	open( OUT, ">>$output") or die;
	print OUT join("\t", ("sample", "5chr_total", "7chr_total",  @chrs)), "\n";
}
foreach my $bam_file(@bam_lists){
	my $bam_input = File::Spec->catfile( $bam_dir, $bam_file);
	die "$bam_file do not exists\n\n" unless (-e $bam_input);
	if ( $bam_file =~ /(\S+)\.bam$/) {
		my $bam_pre = $1;
		my @nums    = ();
		for my $chr( @chrs){
			my $cmd = "samtools view -c $bam_input  $chr";
			if ($debug) {
				print STDERR $cmd, "\n";#code
			}
			
			unless($debug){
				my $n   = `$cmd`;
				chomp $n;
				push @nums, $n;
			}
			
			
		}
		unless($debug){
			my $sum7 = get_sum(@nums);
			my $sum5  = get_sum(@nums[0..4]);
			print OUT join("\t", ($bam_pre, $sum5, $sum7, @nums)), "\n";
		}
	}else{
		die "bam";
	}
}

unless($debug){
	close OUT;
}
exit;

sub get_sum{
	my @tmp = @_;
	my $s = 0;
	foreach my $i(@tmp){
		$s+=$i;
	}
	return $s;
}
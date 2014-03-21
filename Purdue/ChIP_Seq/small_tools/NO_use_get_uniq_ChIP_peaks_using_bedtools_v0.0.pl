#!/usr/bin/perl -w

# Mar 7, 2014
# Now we have ChIP-Seq data that each samples has a Input,
# Compare ChIP and Input can get peaks
# For mutants and WT, we input the two peaks and want to get two
# output files with uniq peaks
# 

use strict;
use File::Spec;

print STDERR "bedtools gives wrong ouput \n\n NO USE\n\n";

my $debug = 0;
if ($debug){
	print STDERR "debug: $debug\n";
}

#print STDERR "input must have a head with Start in it.\n\n";

my $usage = "$0\n <bed_file1_mut> <bed_file2_WT> <outdir> <label>\n\n";
die $usage unless (@ARGV == 4);

my $in1  = shift or die;
my $in2  = shift or die;
#my $allowed_gap = 0;
#my $allowed_gap = shift;

my $outdir = shift or die;
my $label  = shift or die;

die "wrong files\n\n" unless (-e $in1 and -e $in2);
die "no outdir" unless (-d $outdir);

my $tmp1 = File::Spec->catfile($outdir, "tmp.001.gff");
die if(-e $tmp1);
my $tmp2 = File::Spec->catfile($outdir, "tmp.002.gff");
die if(-e $tmp2);


#get_beds($in1, $tmp1);
#get_beds($in2, $tmp2);
get_gff($in1, $tmp1);
get_gff($in2, $tmp2);


my $uniq_tmp1 = File::Spec->catfile($outdir, "uniq_tmp.001.gff");
die if(-e $uniq_tmp1);
my $uniq_tmp2 = File::Spec->catfile($outdir, "uniq_tmp.002.gff");
die if(-e $uniq_tmp2);

#-v	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
#bedtools intersect [OPTIONS] -a <BED/BAM/GFF/VCF> -b <BED/BAM/GFF/VCF>
get_uniq_gff($tmp1, $tmp2, $uniq_tmp1);
get_uniq_gff($tmp2, $tmp1, $uniq_tmp2);

get_final_files($in1, $uniq_tmp1, $outdir, "" );
get_final_files($in2, $uniq_tmp2, $outdir, $label );

my $rm_cmd = "rm -f $tmp1 $tmp2 $uniq_tmp1 $uniq_tmp2";
#`$rm_cmd`;

exit;

#get_beds($in1, $tmp1);
sub get_beds{
	my ($input, $output) = @_;
	die unless (-e $input);
	die if (-e $output);
	open(IN , $input) or die;
	open(OUT , ">$output") or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		next if ($a[0] =~ /^chr$/i);
		print OUT join("\t", @a[0..2]), "\n";
	}
	close IN;
	close OUT;
}


sub get_gff{
	my ($input, $output) = @_;
	die unless (-e $input);
	die if (-e $output);
	open(IN , $input) or die;
	open(OUT , ">$output") or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		next if ($a[0] =~ /^chr$/i);
		#
		#my @dots = (".") x 4;
		print OUT join("\t", $a[0], ".", ".", @a[1..2], ".", "+", ".", "."), "\n";
	}
	close IN;
	close OUT;
}

#get_uniq_beds($in1, $in2, $uniq_tmp1);
sub get_uniq_gff{
	my ($bedA, $bedB, $output) = @_;
	die unless (-e $bedA and -e $bedB);
	my $cmd = "bedtools intersect -v -a $bedA -b $bedB > $output";
	print STDERR $cmd, "\n\n";
	if ( ! $debug) {
		`$cmd`;
	}
}

#get_final_files($in2, $uniq_tmp2, $outdir, $label )
sub get_final_files{
	my ($input, $bed, $outdir_sub, $label_sub) = @_;
	die unless (-e $input and -e $bed and -d $outdir_sub);
	my %records ;
	
	my $output = "";
	my ($volume,$directories,$file) =  File::Spec->splitpath( $input );
	if ($file =~ /(\S+)\.txt/) {
		my $pre = $1;
		
		if ($label_sub eq "" ) {
			$output = File::Spec->catfile($outdir, $pre . "_uniq.txt");
		}else{
			$output = File::Spec->catfile($outdir, $pre . "_uniq_" . $label_sub . ".txt");
		}
		
		
		
		die if(-e $output);
	}
	
	
	open( IN, $bed  ) or die;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $tmp = join("_", @a[0,3,4]);
		$records{$tmp} = 1;
	}
	close IN;
	
	open( IN, $input) or die;
	open(OUT, ">>$output" ) or die;
	my $h = <IN>;
	print OUT $h;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		my $tmp = join("_", @a[0..2]);
		if( defined $records{$tmp} ){
			print OUT $_, "\n";
		}
	}
	
	close IN;
	close OUT;	
}

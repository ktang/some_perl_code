#!/usr/bin/perl -w

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $index_db = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome";

#my $ref_fas = 	"/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa";
#die unless (-e $ref_fas);

#my $ref_fai = 	"/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa.fai";
#die unless(-e $ref_fai);


my %phred_label = ( 33=> "--phred33", 64=>"--phred64");
# --phred33          qualities are Phred+33 (default)
# --phred64          qualities are Phred+64

my $usage = "\n$0 \n\n<indir> <inpre> <outdir> <outpre> <phred_33_64>\n\n";
die $usage unless(@ARGV == 5);

my $indir     = shift or die;
my $inpre     = shift or die;
my $outdir    = shift or die;
my $outpre    = shift or die;
my $phred_num = shift or die;
die unless (defined $phred_label{$phred_num});

my $fq1 = File::Spec->catfile($indir, $inpre . "_R1.fq.gz");
my $fq2 = File::Spec->catfile($indir, $inpre . "_R2.fq.gz");

die "input do not exist\n\n" unless (-e $fq1 and -e $fq2);

my $log_file = File::Spec->catfile($outdir, $outpre . "_log.txt");
die "log file exist\n\n" if(-e $log_file);

##########
# bowtie2
#############
my $bam_file = File::Spec->catfile($outdir, $outpre . "_bowtie2_CP.bam");
die if(-e $bam_file);
my $k = 2;

my $p = 4;
my $cmd_bowtie2 = "bowtie2 -t -p $p -k $k $phred_label{$phred_num} --sensitive  -I 0 -X 1000 -x $index_db -1 $fq1 -2 $fq2 2>> $log_file | perl -lane 'if(/^@/) {print} else{print if (\$F[4] == 255 and \$F[-1] eq \"YT:Z:CP\")}' | samtools view -b -S - > $bam_file";

print STDERR $cmd_bowtie2, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_bowtie2, "\n";
	close OUT;
	`$cmd_bowtie2`;
}

#######
# sort
########
my $sorted_bam_pre = File::Spec->catfile($outdir, $outpre . "_bowtie2_CP_sort");
my $sorted_bam = File::Spec->catfile($outdir, $outpre . "_bowtie2_CP_sort.bam");
if(!$debug){
	die unless(-e $bam_file);
	die if(-e $sorted_bam);
}
my $cmd_sort = "samtools sort $bam_file $sorted_bam_pre";

print STDERR $cmd_sort, "\n\n";

unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_sort, "\n\n\n";
	close OUT;
	`$cmd_sort`;
}

#######
# rmdup
########
my $bam_rmdup = File::Spec->catfile($outdir, $outpre . "_bowtie2_CP_sort_rmdup.bam");

if(!$debug){
	die unless(-e $sorted_bam);
	die if(-e $bam_rmdup);
}
my $cmd_rmdup = "samtools rmdup $sorted_bam $bam_rmdup 2>> $log_file";
print STDERR $cmd_rmdup, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_rmdup, "\n";
	close OUT;
	`$cmd_rmdup`;
}

######
# index
##########
my $cmd_index = "samtools index $bam_rmdup";
print STDERR $cmd_index, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_index, "\n\n";
	close OUT;
	`$cmd_index`;
}

exit;

=head
if(!$debug){
	die unless(-e $);
	die if(-e $);
}

print STDERR $cmd_, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_, "\n\n";
	close OUT;
	`$cmd_`;
}



=cut

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
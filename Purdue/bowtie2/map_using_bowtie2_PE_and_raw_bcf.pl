#!/usr/bin/perl -w

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $post_fix = "XX";

my $index_db = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome";
my $genome_fa_file = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa";
die unless (-e $genome_fa_file);

my %phred_label = ( 33=> "--phred33", 64=>"--phred64");
# --phred33          qualities are Phred+33 (default)  L - Illumina 1.8+ Phred+33
# --phred64          qualities are Phred+64 J - Illumina 1.5+ Phred+64

my $usage = "\n$0 \n\n<indir> <inpre> <outdir> <outpre> <phred_33_64>\n\n";
die $usage unless(@ARGV == 5);

my $indir     = shift or die;
my $inpre     = shift or die;
my $outdir    = shift or die;
my $outpre    = shift or die;
my $phred_num = shift or die;
die unless (defined $phred_label{$phred_num});

#my $fq1 = File::Spec->catfile($indir, $inpre . "_R1.fq.gz");
#my $fq2 = File::Spec->catfile($indir, $inpre . "_R2.fq.gz");
#die "input do not exist\n\n" unless (-e $fq1 and -e $fq2);

my @fqs_R1 = get_fq_files("_R1.fq.gz", $indir, $inpre );
my @fqs_R2 = get_fq_files("_R2.fq.gz", $indir, $inpre );

my $fq1 = join(",", @fqs_R1);
my $fq2 = join(",", @fqs_R2);

opendir(INDIR, $indir) or die;
my @fqs1 = grep /R1.fq.gz/, grep /$inpre/, readdir INDIR;
closedir INDIR;

my $log_file = File::Spec->catfile($outdir, $outpre . "_".$post_fix. "_log.txt");
die "log file exist\n\n" if(-e $log_file);

##########
# bowtie2
#############
my $bam_file = File::Spec->catfile($outdir, $outpre . "_bowtie2_" . $post_fix . ".bam");
die if(-e $bam_file);
#my $k = 10;
my $MAPQ_cutoff = 20;
my $p = 2;
my $trim5 = 8;
my $trim3  = 12;
#my $cmd_bowtie2 = "bowtie2 -t -p $p -k $k $phred_label{$phred_num} --very-sensitive --no-mixed --no-discordant -I 0 -X 1000 -x $index_db -1 $fq1 -2 $fq2 2>> $log_file | perl -lane ' print if( \$F[2] ne \"*\")' | samtools view -b -S - > $bam_file";
#my $cmd_bowtie2 = "bowtie2 -t -p $p -k $k $phred_label{$phred_num} --very-sensitive --no-mixed --no-discordant -I 0 -X 1000 -x $index_db -1 $fq1 -2 $fq2 2>> $log_file |samtools view -b -S - > $bam_file";
#my $cmd_bowtie2 = "bowtie2 -t -p $p -k $k $phred_label{$phred_num} --very-sensitive  -I 0 -X 1000 -x $index_db -1 $fq1 -2 $fq2 2>> $log_file |samtools view -b -S - > $bam_file";

my $cmd_bowtie2 = "bowtie2 -t -p $p $phred_label{$phred_num} --trim5 $trim5 --trim3 $trim3 --very-sensitive  -I 0 -X 1000 -x $index_db -1 $fq1 -2 $fq2 2>> $log_file | perl -lane ' print if( \$F[4] >= $MAPQ_cutoff or /^@/ )' | samtools view -b -S - > $bam_file";

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
my $sorted_bam_pre = File::Spec->catfile($outdir, $outpre . "_bowtie2_" . $post_fix . "_sort"  );
my $sorted_bam     = File::Spec->catfile($outdir,  $outpre . "_bowtie2_" . $post_fix . "_sort.bam" );
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
my $bam_rmdup = File::Spec->catfile($outdir,$outpre . "_bowtie2_" . $post_fix . "_sort_rmdup.bam" );

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
my $cmd_index_raw = "samtools index $sorted_bam";
print STDERR $cmd_index_raw, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_index_raw, "\n\n";
	close OUT;
	`$cmd_index_raw`;
}


my $cmd_index = "samtools index $bam_rmdup";
print STDERR $cmd_index, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_index, "\n\n";
	close OUT;
	`$cmd_index`;
}


######
# bcf
##########
#time samtools mpileup -uf /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa rmdup_bam_ln/4/28_XZ10_nrpe_E_bowtie2_XH_Dec24_sort_rmdup.bam | bcftools view -bvcg - >
#2_rmdup_raw_bcf/28_XZ10_nrpe_E_bowtie2_XH_Dec24_sort_rmdup.raw.bcf

#my $genome_fa_file = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa";

my $raw_bcf_file = File::Spec->catfile($outdir,$outpre . "_bowtie2_" . $post_fix . "_sort_rmdup.raw.bcf" );

if(!$debug){
	die unless(-e $bam_rmdup);
	die if(-e $raw_bcf_file);
}
my $cmd_bcf = " samtools mpileup -uf $genome_fa_file $bam_rmdup |  bcftools view -bvcg - > $raw_bcf_file";
print STDERR $cmd_bcf, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT "\n\n", $cmd_bcf, "\n\n";
	close OUT;
	`$cmd_bcf`;
}

exit;

##################################
#my @fqs_R2 = get_fq_files("_R2.fq.gz", $indir, $inpre );
#opendir(INDIR, $indir) or die;
#my @fqs1 = grep /R1.fq.gz/, grep /$inpre/, readdir INDIR;
#closedir INDIR;
sub get_fq_files{
	my ($num_tag, $dir, $pre) = @_;
	die unless (-d $dir);
	opendir(INDIR, $dir) or die;
	my @names = grep /$num_tag/, grep /$dir/, readdir INDIR;
	closedir INDIR;
	
	my @files = ();
	foreach my $t (@names){
		my $f = File::Spec->catfile($dir,$t);
		die $t unless (-e $f);
		push @files, $f;
	}
	
	return @files;	
}
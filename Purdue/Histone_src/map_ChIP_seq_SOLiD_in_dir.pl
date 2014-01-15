#!/usr/bin/perl -w

use strict;
use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);


my $fastq_dump = "/Users/tang58/Software/NCBI/sratoolkit.2.1.10-mac64/bin/fastq-dump";
die unless (-e $fastq_dump);

my $index_db = "/Users/tang58/DataBase/TAIR_Col0_genome/index/Bowtie/ColorSpace/a_thaliana_c";

my $bowtie   = "/Users/tang58/Software/Bowtie/bowtie1/bowtie-0.12.8/bowtie";
die unless (-e $bowtie);

my $fai_file = "/Users/tang58/DataBase/TAIR_Col0_genome/index/TAIR10_Col0_7Chr_for_bowtie/TAIR10_Col0_7Chr.fa.fai";
die unless (-e $fai_file);
#####################################################################
my $debug = 1;


my $usage = "$0 <dir> <label>";
die $usage unless(@ARGV == 2);

my $dir = shift or die;
die unless (-d $dir);

my $label = shift or die;

opendir ( DIR, $dir ) or die;

my @files = grep /\.sra$/, readdir DIR;

print STDERR "\n\n";

foreach my $sra_in (@files){
	
	print STDERR "handling $sra_in...\n";
	
	my $pre = "";
	if( $sra_in =~ /(\S+)\.sra$/ ){
		$pre = $1;
	}else{
		die "pre $sra_in\n\n";
	}
	
	my $sra_file = File::Spec->catfile($dir, $sra_in);
	
	my $log_file  = File::Spec->catfile($dir, $pre . "_$label" . "_log.txt");
	
	my $fastq_file =  File::Spec->catfile($dir, $pre . ".fastq");
	
	my $fq = $pre . ".fastq";
	
	die $fastq_file  if(-e $fastq_file);
	
	##########
	#	1 extract fastq
	###############
	my $cmd_extract = "$fastq_dump --defline-qual \"+\" -E --defline-seq \"@\\\$ac_\\\$si\" $sra_file >> $log_file  2>&1 ";
	
	print STDERR $cmd_extract, "\n\n";
	
	if(! $debug){
		open(IN, ">>$log_file") or die;
		print IN $cmd_extract, "\n";
		close(IN);
		
		`$cmd_extract`;
		`mv $fq $dir`;
		die $fastq_file unless (-e $fastq_file);
	}
	
	##########
	#	2 bowtie map
	###############
#time /Users/tang58/Software/Bowtie/bowtie1/bowtie-0.12.8/bowtie -S -v 3 -m 1 -C /Users/tang58/DataBase/TAIR_Col0_genome/index/Bowtie/ColorSpace/a_thaliana_c -q SRR037881.fastq 
#| awk ' {if(  $3 != "*" )  {print $0} }' > SRR037881_H3K27me3_bowtie_v3m1p1.sam   2>>  $log_file
	my $v = 3;
	my $m = 1;
	my $p = 4;
	
	
	my $sam_file =  File::Spec->catfile($dir, $pre . "_$label" . "_bowtie_" . "v" . $v . "m" . $m . "p" .  $p. ".sam");
	die $sam_file  if(-e $sam_file);
	my $cmd_bowtie = "$bowtie -S -v $v -m $m -p $p -C $index_db -q $fastq_file  2>> $log_file | awk \' {if(  \$3 != \"*\" )  {print \$0} }\' >  $sam_file ";
	
	print STDERR $cmd_bowtie, "\n\n";
	if(! $debug){
		open(IN, ">>$log_file") or die;
		print IN "\n\n", $cmd_bowtie, "\n";
		close(IN);
		
		`$cmd_bowtie`;
		die $sam_file unless (-e $sam_file);
	}
	
	
	##########
	#	3 sam2bam
	###############
	my $bam_raw_file = File::Spec->catfile($dir, $pre . "_$label" . "_bowtie_" . "v" . $v . "m" . $m . "p" .  $p. ".bam");
	die $bam_raw_file if (-e $bam_raw_file);
	my $cmd_sam2bam = "samtools view -bt $fai_file $sam_file > $bam_raw_file 2>> $log_file";
	
	print STDERR $cmd_sam2bam, "\n\n";
	if(! $debug){
		open(IN, ">>$log_file") or die;
		print IN "\n\n", $cmd_sam2bam, "\n";
		close(IN);
		
		`$cmd_sam2bam`;
		die $bam_raw_file unless (-e $bam_raw_file);
	}
	
	##########
	#	4 sort bam
	###############
	my $bam_sorted_pre  =  File::Spec->catfile($dir, $pre . "_$label" . "_bowtie_" . "v" . $v . "m" . $m . "p" .  $p . "_sorted");
	my $bam_sorted_file =  File::Spec->catfile($dir, $pre . "_$label" . "_bowtie_" . "v" . $v . "m" . $m . "p" .  $p. "_sorted.bam");
	die $bam_sorted_file  if (-e $bam_sorted_file );
	
	my $cmd_sort = "samtools sort $bam_raw_file $bam_sorted_pre 2>> $log_file";
	
	print STDERR $cmd_sort, "\n\n";
	if(! $debug){
		open(IN, ">>$log_file") or die;
		print IN "\n\n", $cmd_sort, "\n";
		close(IN);
		
		`$cmd_sort`;
		die $bam_sorted_file unless (-e $bam_sorted_file);
	}
	
	##########
	#	5 rmdup
	###############
	
	my $bam_redup = File::Spec->catfile($dir, $pre . "_$label" . "_bowtie_" . "v" . $v . "m" . $m . "p" .  $p. "_sorted_rmdup.bam");
	die $bam_redup if (-e $bam_redup);
	
	my $cmd_redup = "samtools rmdup -s $bam_sorted_file $bam_redup 2>> $log_file ";	
	
	print STDERR $cmd_redup, "\n\n\n\n";
	if(! $debug){
		open(IN, ">>$log_file") or die;
		print IN "\n\n", $cmd_redup, "\n";
		close(IN);
		
		`$cmd_redup`;
		die  $bam_redup unless (-e $bam_redup);
	}
}


exit;


=head
	print STDERR $cmd_, "\n\n";
	if(! $debug){
		`$cmd_`;
		die  $ unless (-e $);
	}

=cut
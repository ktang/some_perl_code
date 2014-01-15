#!/usr/bin/perl -w
=head2
Copyright (C) 2010 Kai Tang
version 1 29 Nov 2010

for Single_end fastq

like a .sh script to batch run trim end, trim adaptor,
soap map, convert to bam, get SNP pileup
This is version 1, all use default parameter
=cut

sub trim_end
{
	my $orignal_file = shift;	
	my $out_dir = shift;
	my $pre = shift;
	my $out = $pre."_trimmdend.fastq";
	my $out_stat = $pre."_trimmdend.stat";
	if (( -e "$out_dir/$out") or (-e "$out_dir/$out_stat") )
	{ die  "$out_dir/$out or $out_dir/$out_stat exists!!\n";  }
	else 
	{
		$trim_cmd = "fastq_quality_trimmer -v -t 20 -l 20 -i $orignal_file -o $out_dir/$out > $out_dir/$out.trim_end_stat";
		print STDERR "$trim_cmd\ntrimming...\n\n";
		`$trim_cmd`;
		print STDERR "trimming finished!\n";
	}
}

sub trim_adapter
{
	my $pre = shift;
	my $dir = shift;
	my $adapter = shift;
	die "no adapter" unless (defined $adapter);
	my $in = $pre."_trimmdend.fastq";
	my $out = $pre."_clean.fastq";
	my $stat =$pre."_trimEndAda.stat";
	if ( (-e "$dir/$out") or (-e "$dir/$stat") or !(-e "$dir/$in") ) 
	{
		die "(-e $dir/$out) or (-e $dir/$stat) or !(-e $dir/$in)"	
	}
	
	else 
	{
		$cmd = "fastx_clipper -a $adapter -l 20 -n -v -i $dir/$in -o $dir/$out > $dir/$stat";
		print STDERR "$cmd\ntrimming adapter\n\n";
		`$cmd`;	
		`rm -f $dir/$in`;
		print STDERR "finished!\n\n";
	}
		
}
################################################################
#
#		main 
#
################################################################

use strict;
my $usage = "$0 <input_dir> <out_dir> <adapter> <map_DB|C24 or Col0>";
die $usage unless(@ARGV == 4);

my ($indir,$outdir,$ada,$db) = @ARGV[0..3];
my $indexFile;
my $ref_fas;
my $ref_fai;

print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";

if ($db eq "C24"){ 
	$indexFile = "/Users/kaitang/Desktop/C24_genome/SOAP_index/C24_TAIR9_5Chr.fas.index";
	$ref_fas = 	"/Users/kaitang/Desktop/C24_genome/SOAP_index/C24_TAIR9_5Chr.fas";
	$ref_fai = "/Users/kaitang/Desktop/C24_genome/SOAP_index/C24_TAIR9_5Chr.fas.fai";
}
elsif ($db eq "Col0") {
	$indexFile = "/Users/kaitang/Desktop/TAIR/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.index";
	$ref_fas = 	"/Users/kaitang/Desktop/TAIR/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas";
	$ref_fai = 	"/Users/kaitang/Desktop/TAIR/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.fai";
}
else {
	die "$usage\n,DB must be C24 or Col0\n"
}

opendir (INDIR, $indir)
  or die "Cannot open dir $indir: $!";
  
 my @files = grep {/\.fastq$/} readdir INDIR;
 
 my $i = 0;
 
 foreach my $file (@files){
 	
 	$i++;
 	my $pre = "pre".$i;
 	
 	if ($file =~ /(\S+)\.fastq$/)
	{$pre = $1;}
	
	# trim end
	trim_end ("$indir/$file",$indir,$pre);
	
	
	#trim adaptor
	trim_adapter($pre, $indir, $ada);
	
	#soap
	my $clean_fastq = $pre."_clean.fastq";
	my $soap_out =  $pre."_SE_n1r0v1M4.soap";
	if (-e "$outdir/$soap_out" or !(-e "$indir/$clean_fastq")){
		print STDERR "$outdir/$soap_out exists!! or $indir/$clean_fastq do not exists!! \n";
		die "soap";
	}

	else{
		my $cmd_SE = "soap -a $indir/$clean_fastq -D $indexFile -o $outdir/$soap_out -n 1 -r 0 -v 1 -M 4";
		print STDERR "$cmd_SE\nsoaping for SE\n\n";
		`$cmd_SE`;
		print STDERR "finished!\n\n";
	}

	#soap2sam	
	my $sam = $pre."_SE_n1r0v1M4.sam";
 
 	if ( (-e "$outdir/$sam") ){
	 	die "die at soap to sam\n$sam\n"	
 	}
 	
 	else{
 		my $soap2sam_cmd1 = "soap2sam.pl $outdir/$soap_out > $outdir/$sam";
 		print STDERR "$soap2sam_cmd1\nsoap2sam converting...\n\n";
 		`$soap2sam_cmd1`;
 		print STDERR "finishes!\n";
 	
 		`rm -f $outdir/$soap_out`;
	 }
	 
	#sam2bam
	my $bam = $pre."_SE_n1r0v1M4.bam";
 
	 if ( (-e "$outdir/$bam") )
 	{
 		die "die at sam2bam\n$bam\n"	
 	}
 	else{
 		my $sam2bam_cmd1 = "samtools view -bt $ref_fai $outdir/$sam > $outdir/$bam";
	 	print  STDERR "$sam2bam_cmd1\nsam to bam converting...\n\n";
 		`$sam2bam_cmd1`;
 		print STDERR "finished!\n\n";
 		`rm -f $outdir/$sam`;
 	}
	
}

# merge bam

opendir (OUTDIR,$outdir)
	or die "cannot open outdir $outdir:$!";

my $merged_bam = "all_merged.bam";
my $merge_bam_cmd = "samtools merge $outdir/$merged_bam";

my @bams = grep /.bam$/ , readdir OUTDIR;

foreach my $bam (@bams){
	$merge_bam_cmd = $merge_bam_cmd." $outdir/$bam";
}

print STDERR "$merge_bam_cmd\nbam merging...\n\n";
`$merge_bam_cmd`;
print STDERR "finished merge!\n\n";

foreach my $bam (@bams){
	`rm -f $outdir/$bam`;
}

##############
#sort bam
##############


 my  $sorted_bam_pre = "all_merged_sorted";
 my $sorted_bam = "all_merged_sorted.bam";
 if (!(-e "$outdir/$sorted_bam") ){
 	my $sort_cmd = "samtools sort $outdir/$merged_bam $outdir/$sorted_bam_pre";
	print STDERR "$sort_cmd\nbam sorting...\n\n";
	`$sort_cmd`;
 	print STDERR "sorting finishes!\n\n";
 	`rm -f $outdir/$merged_bam`;
 }
 
 else{
 	die "sorted bam exists!!\n";
 }
 
 ##############
 #  index
 ##############
 if (-e "$outdir/$sorted_bam"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam";
	 print STDERR "$index_cmd\nindexing\n\n";
	 `$index_cmd`;
	 print STDERR "finished index!\n\n";
 }
 
 else{
 	die "no $outdir/$sorted_bam";
 }
	
	
 # pileup
 my $pileup_file ="SNP.pileup";
 if (!(-e "$outdir/$pileup_file")){
 	my $pileup_cmd = "samtools pileup -f $ref_fas $outdir/$sorted_bam | pileup_filter.pl > $outdir/$pileup_file";
 	print STDERR "$pileup_cmd\npileupping...\n\n\n";
	`$pileup_cmd`;
	print STDERR "finished all!!\n";
 }
 else{
 	die "pileup exists!\n";
 }
exit;
#!/usr/bin/perl -w
=head2

v1.2 partly
already get adaptor trimmed file, do the following work

v1.1(Feb 15, 2011)
get pair use fasta as output, otherwise the memory in not enough.

Copyright (C) 2010 Kai Tang
version 1 27 Nov 2010

like a .sh script to batch run trim end, trim adaptor,
soap map, convert to bam, get SNP pileup
This is version 1, all use default parameter
=cut

################################################################
#		main 
################################################################

print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";

use strict;
my $usage = "$0 <pre> <map_DB|C24 or Col0>";
print STDERR "working directory should be the same directory\n\n";
die $usage unless(@ARGV == 2);

my ($pre, $db) = @ARGV[0..1];

#my ($indir,$outdir,$ada1,$ada2,$db) = @ARGV[0..4];
my $indexFile = "";
my $ref_fas = "";
my $ref_fai = "";

if ($db eq "C24"){ 
	$indexFile = "/Users/tang58/DataBase/C24/index/SOAP/C24_TAIR9_5Chr.fas.index";
	$ref_fas = 	"/Users/tang58/DataBase/C24/C24_TAIR9_5Chr.fas";
	$ref_fai = "/Users/tang58/DataBase/C24/index/SOAP/C24_TAIR9_5Chr.fas.fai";
}
elsif ($db eq "Col0") {
	$indexFile = "/Users/tang58/DataBase/TAIR_Col0_genome/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.index";
	$ref_fas = 	"/Users/tang58/DataBase/TAIR_Col0_genome/5Chr_only_TAIR9_Col0.fas";
	$ref_fai = 	"/Users/tang58/DataBase/TAIR_Col0_genome/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.fai";
}
else {
	die "$usage\n,DB must be C24 or Col0\n"
}

#get pairs
my $clean_pair1 = $pre."_1_pairs.fastq";
my $clean_pair2 = $pre."_2_pairs.fastq";

#soap_pair
my $pe_out = $pre."_PE_m0x1000n1r0v2M4.soap";
my $left_out =  $pre. "_PEleft_m0x1000n1r0v2M4.soap";
if ( (-e $pe_out) or (-e $left_out)){ 
	print STDERR "die from soap_pe\n";  
	die 
}
else{
	#v 1-> v4
	my $soap_pe = "time soap -a $clean_pair1 -b $clean_pair2 -D $indexFile -o $pe_out -2 $left_out -m 0 -x 1000 -n 1 -r 0 -v 2 -M 4";
	print STDERR "$soap_pe\n soaping for pair end...\n";
	`$soap_pe`;
}

#soap_single_end
my $clean_single1 = $pre."_1_single.fastq";
my $clean_single2 = $pre."_2_single.fastq";

my $single_out1 =  $pre."_1_single_n1r0v2M4.soap";
my $single_out2 =  $pre."_2_single_n1r0v2M4.soap";

if (-e $single_out1 or -e $single_out2 ){
	die "sigle_end";
}
else{
	my $cmd_SE1 = "time soap -a $clean_single1 -D $indexFile -o $single_out1 -n 1 -r 0 -v 2 -M 4";
	my $cmd_SE2 = "time soap -a $clean_single2 -D $indexFile -o $single_out2 -n 1 -r 0 -v 2 -M 4";
	print STDERR "$cmd_SE1\n";
	`$cmd_SE1`;
	print STDERR "$cmd_SE2\n";
	`$cmd_SE2`;
}

# soap2sam
 my $PE_sam		= $pre."_PE_m0x1000n1r0v2M4.sam";
 my $left_sam 	= $pre. "_PEleft_m0x1000n1r0v2M4.sam";
 my $SE_sam1 	= $pre."_1_single_n1r0v2M4.sam";
 my $SE_sam2 	= $pre."_2_single_n1r0v2M4.sam";
 
 if (-e $PE_sam  or -e $left_sam  or -e $SE_sam1  or -e $SE_sam2) {
 	die "die at soap to sam\n$PE_sam\n"	
 }
 else{
 	my $soap2sam_cmd1 = "time soap2sam.pl -p $pe_out > $PE_sam";
 	print STDERR "$soap2sam_cmd1\n";
 	`$soap2sam_cmd1`;
 	
 	my $soap2sam_cmd2 = "time soap2sam.pl $left_out > $left_sam";
 	print STDERR "$soap2sam_cmd2\n";
 	`$soap2sam_cmd2`;
 	
	my $soap2sam_cmd3 = "time soap2sam.pl $single_out1 > $SE_sam1";
 	print STDERR "$soap2sam_cmd3\n";
 	`$soap2sam_cmd3`;
 	
 	my $soap2sam_cmd4 = "time soap2sam.pl $single_out2 > $SE_sam2";
 	print STDERR "$soap2sam_cmd4\n";
 	`$soap2sam_cmd4`;
 	
 #	`rm -f $outdir/$pe_out`;
 #	`rm -f $outdir/$left_out`;
 #	`rm -f $outdir/$single_out`;
}
 
 
#sam to bam
my $PE_bam		= $pre."_PE_m0x1000n1r0v2M4.bam";
my $left_bam 	= $pre. "_PEleft_m0x1000n1r0v2M4.bam";
my $SE_bam1 	= $pre."_1_single_n1r0v2M4.bam";
my $SE_bam2 	= $pre."_2_single_n1r0v2M4.bam";
 
 if ( -e $PE_bam ){
 	die "die at sam2bam\n"	
 }
 else{
 	my $sam2bam_cmd1 = "time samtools view -bt $ref_fai $PE_sam > $PE_bam ";
 	print STDERR "$sam2bam_cmd1\n";
 	`$sam2bam_cmd1`;
 #	`rm -f $outdir/$PE_sam`;
 }
 
 if (  -e $left_bam){
 		die "die at sam2bam\n$left_bam\n"	
 }
 else{	
 	my $sam2bam_cmd2 = "time samtools view -bt $ref_fai $left_sam > $left_bam ";
 	print STDERR "$sam2bam_cmd2\n";
 	`$sam2bam_cmd2`;
 	
 #	`rm -f $outdir/$left_sam`;
 }
 
  if ( -e $SE_bam1){
 	die "die at sam2bam\n"	
  }
  else{	
 	my $sam3bam_cmd3 = "time samtools view -bt $ref_fai $SE_sam1 > $SE_bam1 ";
 	print STDERR "$sam3bam_cmd3\n";
 	`$sam3bam_cmd3`;
 	
 #	`rm -f $outdir/$SE_sam`;
 }

if ( -e $SE_bam2){
 	die "die at sam2bam\n"	
 }
 else{	
 	my $sam3bam_cmd4 = "time samtools view -bt $ref_fai $SE_sam2 > $SE_bam2 ";
 	print STDERR "$sam3bam_cmd4\n";
 	`$sam3bam_cmd4`;
 	
 #	`rm -f $outdir/$SE_sam`;
 }
  
# merge bam
my $merged_bam = $pre. "_all_merged.bam";
my $merge_bam_cmd = "time samtools merge $merged_bam $PE_bam $left_bam $SE_bam1 $SE_bam2";
print STDERR "$merge_bam_cmd\n";
`$merge_bam_cmd`;
#`rm -f $outdir/$PE_bam $outdir/$left_bam $outdir/$SE_bam`;

# sort bam
my $sorted_bam_pre = $pre."_all_merged_sorted";
my $sorted_bam = $pre."_all_merged_sorted.bam";
my $sort_cmd = "time samtools sort $merged_bam $sorted_bam_pre ";
print STDERR "$sort_cmd\n";
 `$sort_cmd`;
# `rm -f $outdir/$merged_bam`;
 
# index bam
if (-e $sorted_bam){
	 my $index_cmd = "time samtools index $sorted_bam ";
	 print STDERR "$index_cmd";
	 `$index_cmd`;
}
else{
 	die "no $sorted_bam";
}
 
# pileup
# my $pileup_file = $pre1."_".$pre2.".pileup";
# my $pileup_cmd = "samtools pileup -f $ref_fas $outdir/$sorted_bam | pileup_filter.pl > $outdir/$pileup_file &";
# print STDERR "$pileup_cmd\n pileupping...\n";
# `$pileup_cmd`;
 
sub_end_program();

############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}
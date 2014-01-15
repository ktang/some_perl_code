#!/usr/bin/perl -w
=head2
Copyright (C) 2010 Kai Tang
version 1 27 Nov 2010

like a .sh script to batch run trim end, trim adaptor,
soap map, convert to bam, get SNP pileup
This is version 1, all use default parameter
=cut

#Feb 7,2011 use perl script to trim adaptor.

################################################################
#
#		main 
#
################################################################

print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";



use strict;
my $usage = "$0 <C24|Col0> <dir> <input_soap_output_file>";
die $usage unless(@ARGV == 3);

my ($db, $dir, $ip_file) = @ARGV[0..2]; #, $ctrl_file, $output) = @ARGV[0..4];
my $outdir = $dir;

if ( !(-e "$dir/$ip_file") ){
	die "input wrong!!\n";
}
 
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




# soap2sam
 my $ip_sam = $ip_file.".sam";
 #my $ctrl_sam = $ctrl_file.".sam";
 #my $output_sam = $output.".sam";
 
 if ( -e "$dir/$ip_sam" ){
  #	print STDERR "why	hy\n\n\n";
 	die "die at soap to sam\n$ip_sam\n"	
 }
 else{
 	my $soap2sam_cmd1 = "soap2sam.pl $dir/$ip_file > $dir/$ip_sam";
 	print STDERR "$soap2sam_cmd1\n converting\n";
 	`$soap2sam_cmd1`;
 	
 #	`rm -f $outdir/$pe_out`;
 }
 
 

 
 #sam to bam
 my $ip_bam = $ip_file.".bam";
 #my $ctrl_bam = $ctrl_file.".bam";
 #my $output_bam = $output.".bam";

 #print STDERR "$outdir/$ip_bam", "\n\n\n\n";
 
 if ( -e "$outdir/$ip_bam") {
 	print STDERR "bam exists!!!\n";
 	die "die at sam2bam\n$ip_bam\n";	
 }
 else{
 	my $sam2bam_cmd1 = "samtools view -bt $ref_fai $outdir/$ip_sam > $outdir/$ip_bam";
 	print STDERR "$sam2bam_cmd1\n converting\n";
 	`$sam2bam_cmd1`;
 	
 	`rm -f $outdir/$ip_sam`;
 }
 
 
 
# merge bam
#my $merged_bam = $pre1."_".$pre2."_all_merged.bam";
#my $merge_bam_cmd = "samtools merge $outdir/$merged_bam $outdir/$PE_bam $outdir/$left_bam $outdir/$SE_bam";
#print STDERR "$merge_bam_cmd\n merging bam\n";
#`$merge_bam_cmd`;
# `rm -f $outdir/$PE_bam $outdir/$left_bam $outdir/$SE_bam`;

# sort bam ip
 my $sorted_bam_pre_ip = $ip_file."_sorted";
 my $sorted_bam_ip = $ip_file."_sorted.bam";
 
 if (-e "$dir/$sorted_bam_ip"){
 	die "sorted_bam exists!!!\n\n";
 }else{
 	my $sort_cmd_ip = "samtools sort $outdir/$ip_bam $outdir/$sorted_bam_pre_ip";
 	print STDERR "$sort_cmd_ip\n sorting...\n";
 	`$sort_cmd_ip`;
 	`rm -f $outdir/$ip_bam`;
 }
 
 # index bam ip
 
 if (-e "$outdir/$sorted_bam_ip"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam_ip";
	 print STDERR "$index_cmd\nindexing\n";
	 `$index_cmd`;
 }
 
 else{
 	die "no $outdir/$sorted_bam_ip";
 }
  
  
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


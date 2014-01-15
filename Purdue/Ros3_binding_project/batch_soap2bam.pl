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
my $usage = "$0 <C24|Col0> <dir> <IP_file_name> <ctrl_file_name> <output_uniq_name>";
die $usage unless(@ARGV == 5);

my ($db, $dir, $ip_file, $ctrl_file, $output) = @ARGV[0..4];
my $outdir = $dir;

if ( !(-e "$dir/$ip_file") or !(-e "$dir/$ctrl_file") or (-e "$dir/$output") ){
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



my $cmd_uniq = "time ./get_Soap_map_exclude_ctrl_v1.2.pl -c $dir/$ctrl_file -t $dir/$ip_file -o $dir/$output";
if (-e "$dir/$output") {die "$dir/$output exists!!!"}

print STDERR $cmd_uniq, "\n";
system("$cmd_uniq");

# soap2sam
 my $ip_sam = $ip_file.".sam";
 my $ctrl_sam = $ctrl_file.".sam";
 my $output_sam = $output.".sam";
 
 if ( (-e "$dir/$ip_sam") )
 {
 	die "die at soap to sam\n$ip_sam\n"	
 }
 else{
 	my $soap2sam_cmd1 = "soap2sam.pl $dir/$ip_file > $dir/$ip_sam";
 	print STDERR "$soap2sam_cmd1\n converting\n";
 	`$soap2sam_cmd1`;
 	
 #	`rm -f $outdir/$pe_out`;
 }
 
  if (  (-e "$dir/$ctrl_sam") )
 {
 	die "die at soap to sam\n$ctrl_sam\n"	
 }
 else{	
 	my $soap2sam_cmd2 = "soap2sam.pl $dir/$ctrl_file > $dir/$ctrl_sam";
 	print STDERR "$soap2sam_cmd2\n converting\n";
 	`$soap2sam_cmd2`;
 	
 #	`rm -f $outdir/$left_out`;
 }
 
  if ( (-e "$dir/$output_sam"))
 {
 	die "die at soap to sam\n$output_sam\n"	
 }
 else{	
 	my $soap2sam_cmd3 = "soap2sam.pl $dir/$output > $dir/$output_sam";
 	print STDERR "$soap2sam_cmd3\n converting\n";
 	`$soap2sam_cmd3`;
 	
 #	`rm -f $outdir/$single_out`;
 }
 

 
 #sam to bam
 my $ip_bam = $ip_file.".bam";
 my $ctrl_bam = $ctrl_file.".bam";
 my $output_bam = $output.".bam";

 
 if ( (-e "$outdir/$ip_bam") )
 {
 	die "die at sam2bam\n$ip_bam\n"	
 }
 else{
 	my $sam2bam_cmd1 = "samtools view -bt $ref_fai $outdir/$ip_sam > $outdir/$ip_bam";
 	print STDERR "$sam2bam_cmd1\n converting\n";
 	`$sam2bam_cmd1`;
 	
 	`rm -f $outdir/$ip_sam`;
 }
 
  if (  (-e "$outdir/$ctrl_bam") )
 {
 	die "die at sam2bam\n$ctrl_bam\n"	
 }
 else{	
 	my $sam2bam_cmd2 = "samtools view -bt $ref_fai $outdir/$ctrl_sam > $outdir/$ctrl_bam";
 	print STDERR "$sam2bam_cmd2\n converting\n";
 	`$sam2bam_cmd2`;
 	
 	`rm -f $outdir/$ctrl_sam`;
 }
 
  if ( (-e "$outdir/$output_bam"))
 {
 	die "die at sam2bam\n$output_bam\n"	
 }
 else{	
 	my $sam3bam_cmd3 = "samtools view -bt $ref_fai $outdir/$output_sam > $outdir/$output_bam";
 	print STDERR "$sam3bam_cmd3\n converting\n";
 	`$sam3bam_cmd3`;
 	
 	`rm -f $outdir/$output_sam`;
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
 my $sort_cmd_ip = "samtools sort $outdir/$ip_bam $outdir/$sorted_bam_pre_ip";
 print STDERR "$sort_cmd_ip\n sorting...\n";
 `$sort_cmd_ip`;
 `rm -f $outdir/$ip_bam`;

# sort bam ctrl
my $sorted_bam_pre_ctrl = $ctrl_file."_sorted";
 my $sorted_bam_ctrl = $ctrl_file."_sorted.bam";
 my $sort_cmd_ctrl = "samtools sort $outdir/$ctrl_bam $outdir/$sorted_bam_pre_ctrl";
 print STDERR "$sort_cmd_ctrl\n sorting...\n";
 `$sort_cmd_ctrl`;
 `rm -f $outdir/$ctrl_bam`;

# sort bam output
my $sorted_bam_pre_out = $output."_sorted";
 my $sorted_bam_out = $output."_sorted.bam";
 my $sort_cmd_out = "samtools sort $outdir/$output_bam $outdir/$sorted_bam_pre_out";
 print STDERR "$sort_cmd_out\n sorting...\n";
 `$sort_cmd_out`;
 `rm -f $outdir/$output_bam`;
 
 
 # index bam ip
 
 if (-e "$outdir/$sorted_bam_ip"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam_ip";
	 print STDERR "$index_cmd\nindexing\n";
	 `$index_cmd`;
 }
 
 else{
 	die "no $outdir/$sorted_bam_ip";
 }
 # index bam ctrl
 
 if (-e "$outdir/$sorted_bam_ctrl"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam_ctrl";
	 print STDERR "$index_cmd\nindexing\n";
	 `$index_cmd`;
 }
 
 else{
 	die "no $outdir/$sorted_bam_ctrl";
 }
 # index bam out
 
 if (-e "$outdir/$sorted_bam_out"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam_out";
	 print STDERR "$index_cmd\nindexing\n";
	 `$index_cmd`;
 }
 
 else{
 	die "no $outdir/$sorted_bam_out";
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


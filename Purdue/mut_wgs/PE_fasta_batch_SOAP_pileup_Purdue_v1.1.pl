#!/usr/bin/perl -w
=head2
v1.1(Feb 15, 2011)
get pair use fasta as output, otherwise the memory in not enough.

Copyright (C) 2010 Kai Tang
version 1 27 Nov 2010

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
		print STDERR "$trim_cmd\n\n trimming\n\n";
		`$trim_cmd`;
	}
}

=head
sub trim_adapter
{
	my $pre = shift;
	my $dir = shift;
	my $adapter = shift;
	die "no adapter" unless (defined $adapter);
	my $in = $pre."_trimmdend.fastq";
	my $out = $pre."_trimEndAda.fastq";
	my $stat =$pre."_trimEndAda.stat";
	if ( (-e "$dir/$out") or (-e "$dir/$stat") or !(-e "$dir/$in") ) 
	{
		die "(-e $dir/$out) or (-e $dir/$stat) or !(-e $dir/$in)"	
	}
	
	else 
	{
		$cmd = "fastx_clipper -a $adapter -l 20 -n -v -i $dir/$in -o $dir/$out > $dir/$stat";
		print STDERR "$cmd \n\n trimming adapter\n\n";
		`$cmd`;	
		`rm -f $dir/$in`;
	}
		
}
=cut
#Feb 7,2011 use perl script to trim adaptor.
sub trim_adapter
{
	my $pre = shift;
	my $dir = shift;
	my $adapter = shift;
	die "no adapter" unless (defined $adapter);
	my $in = $pre."_trimmdend.fastq";
	my $out = $pre."_trimEndAda.fastq";
	my $stat =$pre."_trimEndAda.stat";
	if ( (-e "$dir/$out") or (-e "$dir/$stat") or !(-e "$dir/$in") ) 
	{
		die "(-e $dir/$out) or (-e $dir/$stat) or !(-e $dir/$in)"	
	}
	
	else 
	{
#		$cmd = "fastx_clipper -a $adapter -l 20 -n -v -i $dir/$in -o $dir/$out > $dir/$stat";
		$cmd = "/Users/tang58/scripts_all/perl_code/Purdue/trim_adaptor_Kai_v1.1.pl --adaptor $adapter --input $dir/$in  --output $dir/$out --minLen 15 > $dir/$stat";

		print STDERR "$cmd \n\n trimming adapter\n\n";
		`$cmd`;	
		`rm -f $dir/$in`;
	}
		
}
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
my $usage = "$0 <input_dir> <out_dir> <adapter1> <adapter2> <map_DB|C24 or Col0>";
die $usage unless(@ARGV == 5);

my ($indir,$outdir,$ada1,$ada2,$db) = @ARGV[0..4];
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

opendir (INDIR, $indir)
  or die "Cannot open dir $indir: $!";
  
 my @files = grep {/\.fastq$/} readdir INDIR;
 
 if ($#files != 1)
	{print STDERR join "\t",readdir INDIR;
		die "not exact two files!";}

my ($pre1,$pre2) = ("pre1","pre2");	
	
if ($files[0] =~ /(\S+)\.fastq$/)
{$pre1 = $1;}

if ($files[1] =~ /(\S+)\.fastq$/)
{$pre2 = $1;}



# trim end
trim_end ("$indir/$files[0]",$indir,$pre1);
trim_end ("$indir/$files[1]",$indir,$pre2);

#trim adaptor

trim_adapter ($pre1, $indir, $ada1);
trim_adapter ($pre2, $indir, $ada2);

#get pairs
my $in1 = $pre1."_trimEndAda.fastq";
my $in2 = $pre2."_trimEndAda.fastq";

my $clean_pair1 = $pre1."_clean.fasta";
my $clean_pair2 = $pre2."_clean.fasta";
my $single_fastq  = $pre1."_".$pre2."_cleanSingle.fasta";
my $pair_stat = $pre1."_".$pre2.".stat";

if ( !(-e "$indir/$in1") or !(-e "$indir/$in2") or (-e "$indir/$clean_pair1") or 
		(-e "$indir/$clean_pair2") or (-e "$indir/$single_fastq") or  (-e "$indir/$pair_stat") )
{die "!(-e $indir/$in1) or !(-e $indir/$in2) or (-e $indir/$clean_pair1) or (-e $indir/$clean_pair2) or (-e $indir/$single_fastq) or  (-e $indir/$pair_stat) "}

else 
{
	my $cmd_pair = "time /Users/tang58/scripts_all/perl_code/Purdue/mut_wgs/find_pairs_for_pair_end_reads_Purdue_v1.5.pl $indir/$in1 $indir/$in2 $indir/$clean_pair1 $indir/$clean_pair2 $indir/$single_fastq > $indir/$pair_stat";
	
	print STDERR "$cmd_pair\n\n getting pair\n";
	`$cmd_pair`; 
#	`rm -f $indir/$in1 $indir/$in2`;
}

#soap_pair
my $pe_out = $pre1."_".$pre2."_PE_m0X1000n1r0v4M4.soap";
my $left_out =  $pre1."_".$pre2."_PEleft_m0X1000n1r0v4M4.soap";
if ( (-e "$outdir/$pe_out") or (-e "$outdir/$left_out"))
{ print STDERR "die from soap_pe\n";  die ;}
else
{
	#v 1-> v4
	my $soap_pe = "soap -a $indir/$clean_pair1 -b $indir/$clean_pair2 -D $indexFile -o $outdir/$pe_out -2 $outdir/$left_out -m 0 -x 1000 -n 1 -r 0 -v 4 -M 4";
	print STDERR "$soap_pe\n soaping for pair end...\n";
	`$soap_pe`;
}
#soap_single_end
my $single_out =  $pre1."_".$pre2."_SE_n1r0v4M4.soap";
if (-e "$outdir/$single_out"){
	die "sigle_end";
}
else{
	my $cmd_SE = "soap -a $indir/$single_fastq -D $indexFile -o $outdir/$single_out -n 1 -r 0 -v 4 -M 4";
	print STDERR "$cmd_SE\n soaping for SE\n";
	`$cmd_SE`;
}
# soap2sam
 my $PE_sam = $pre1."_".$pre2."_PE_m0X1000n1r0v4M4.sam";
 my $left_sam = $pre1."_".$pre2."_PEleft_m0X1000n1r0v4M4.sam";
 my $SE_sam = $pre1."_".$pre2."_SE_n1r0v4M4.sam";
 
 if ( (-e "$outdir/$PE_sam") )
 {
 	die "die at soap to sam\n$PE_sam\n"	
 }
 else{
 	my $soap2sam_cmd1 = "soap2sam.pl -p $outdir/$pe_out > $outdir/$PE_sam";
 	print STDERR "$soap2sam_cmd1\n converting\n";
 	`$soap2sam_cmd1`;
 	
 #	`rm -f $outdir/$pe_out`;
 }
 
  if (  (-e "$outdir/$left_sam") )
 {
 	die "die at soap to sam\n$left_sam\n"	
 }
 else{	
 	my $soap2sam_cmd2 = "soap2sam.pl $outdir/$left_out > $outdir/$left_sam";
 	print STDERR "$soap2sam_cmd2\n converting\n";
 	`$soap2sam_cmd2`;
 	
 #	`rm -f $outdir/$left_out`;
 }
 
  if ( (-e "$outdir/$SE_sam"))
 {
 	die "die at soap to sam\n$SE_sam\n"	
 }
 else{	
 	my $soap2sam_cmd3 = "soap2sam.pl $outdir/$single_out > $outdir/$SE_sam";
 	print STDERR "$soap2sam_cmd3\n converting\n";
 	`$soap2sam_cmd3`;
 	
 #	`rm -f $outdir/$single_out`;
 }
 

 
 #sam to bam
 my $PE_bam = $pre1."_".$pre2."_PE_m0X1000n1r0v4M4.bam";
 my $left_bam = $pre1."_".$pre2."_PEleft_m0X1000n1r0v4M4.bam";
 my $SE_bam = $pre1."_".$pre2."_SE_n1r0v4M4.bam";
 
 if ( (-e "$outdir/$PE_bam") )
 {
 	die "die at sam2bam\n$PE_bam\n"	
 }
 else{
 	my $sam2bam_cmd1 = "samtools view -bt $ref_fai $outdir/$PE_sam > $outdir/$PE_bam";
 	print STDERR "$sam2bam_cmd1\n converting\n";
 	`$sam2bam_cmd1`;
 	
 	`rm -f $outdir/$PE_sam`;
 }
 
  if (  (-e "$outdir/$left_bam") )
 {
 	die "die at sam2bam\n$left_bam\n"	
 }
 else{	
 	my $sam2bam_cmd2 = "samtools view -bt $ref_fai $outdir/$left_sam > $outdir/$left_bam";
 	print STDERR "$sam2bam_cmd2\n converting\n";
 	`$sam2bam_cmd2`;
 	
 	`rm -f $outdir/$left_sam`;
 }
 
  if ( (-e "$outdir/$SE_bam"))
 {
 	die "die at sam2bam\n$SE_bam\n"	
 }
 else{	
 	my $sam3bam_cmd3 = "samtools view -bt $ref_fai $outdir/$SE_sam > $outdir/$SE_bam";
 	print STDERR "$sam3bam_cmd3\n converting\n";
 	`$sam3bam_cmd3`;
 	
 	`rm -f $outdir/$SE_sam`;
 }
 

 
# merge bam
my $merged_bam = $pre1."_".$pre2."_all_merged.bam";
my $merge_bam_cmd = "samtools merge $outdir/$merged_bam $outdir/$PE_bam $outdir/$left_bam $outdir/$SE_bam";
print STDERR "$merge_bam_cmd\n merging bam\n";
`$merge_bam_cmd`;
 `rm -f $outdir/$PE_bam $outdir/$left_bam $outdir/$SE_bam`;

# sort bam
my $sorted_bam_pre = $pre1."_".$pre2."_all_merged_sorted";
 my $sorted_bam = $pre1."_".$pre2."_all_merged_sorted.bam";
 my $sort_cmd = "samtools sort $outdir/$merged_bam $outdir/$sorted_bam_pre";
 print STDERR "$sort_cmd\n sorting...\n";
 `$sort_cmd`;
 `rm -f $outdir/$merged_bam`;
 
 
 # index bam
 
 if (-e "$outdir/$sorted_bam"){
	 my $index_cmd = "samtools index $outdir/$sorted_bam";
	 print STDERR "$index_cmd\nindexing\n";
	 `$index_cmd`;
 }
 
 else{
 	die "no $outdir/$sorted_bam";
 }
 
 # pileup
 my $pileup_file = $pre1."_".$pre2.".pileup";
 my $pileup_cmd = "samtools pileup -f $ref_fas $outdir/$sorted_bam | pileup_filter.pl > $outdir/$pileup_file";
 print STDERR "$pileup_cmd\n pileupping...\n";
 `$pileup_cmd`;
 
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


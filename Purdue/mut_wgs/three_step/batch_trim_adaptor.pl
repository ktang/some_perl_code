#!/usr/bin/perl -w
=head2

part of PE_fasta_batch_SOAP_pileup_Purdue_v1.1.pl
################################################
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
my $usage = "$0 <input_dir> <adapter1> <adapter2>";
die $usage unless(@ARGV == 3);

my ($indir,$ada1,$ada2) = @ARGV[0..2];
my $indexFile = "";
my $ref_fas = "";
my $ref_fai = "";


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
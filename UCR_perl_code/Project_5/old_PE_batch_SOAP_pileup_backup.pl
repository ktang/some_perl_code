#!/usr/bin/perl -w
=head2
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
	if (( -e $out_dir/$out) or (-e $out_dir/$out_stat) )
	{ die  "$out_dir/$out or $out_dir/$out_stat exists!!\n";  }
	else 
	{
		$trim_cmd = "fastq_quality_trimmer -v -t 20 -l 20 -i $orignal_file -o $out_dir/$out > $out.trim_end_stat";
		print "$trim_cmd\n\n trimming\n\n";
		`$trim_cmd`;
	}
}

sub trim_adapter
{
	my $pre = shift;
	my $dir = shift;
	my $adapter = shift;
	die "no adapter" unless (defined $adapter);
	my $in = $pre."_trimmdend.fastq";
	my $out = $pre."_trimEndAda.fastq";
	my $stat =$pre."_trimEndAda.stat";
	if ( (-e $dir/$out) or (-e $dir/$stat) or !(-e $dir/$in) ) 
	{
		die "(-e $dir/$out) or (-e $dir/$stat) or !(-e $dir/$in)"	
	}
	
	else 
	{
		$cmd = "fastx_clipper -a $adapter -l 20 -n -v -i $dir/$in -o $dir/$out > $dir/$stat";
		print "$cmd \n\n trimming adapter\n\n";
		`$cmd`;	
	}
		
}
################################################################
#
#		main 
#
################################################################

use strict;
my $usage = "$0 <input_dir> <out_dir> <adapter1> <adapter2> <map_DB|C24 or Col0>";
die $usage unless(@ARGV == 5);

my ($indir,$outdir,$ada1,$ada2,$db) = @ARGV[0..4];
my $indexFile;

if ($db eq "C24"){ $indexFile = "/Users/kaitang/Desktop/C24_genome/SOAP_index/C24_TAIR9_5Chr.fas.index"}
elsif ($db eq "Col0") {$indexFile = "/Users/kaitang/Desktop/TAIR/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.index"}
else {die "$usage\n,DB must be C24 or Col0\n"}

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
trim_end ($indir/$files[0],$indir,$pre1);
trim_end ($indir/$files[1],$indir,$pre2);

#trim adaptor

trim_adapter ($pre1, $indir, $ada1);
trim_adapter ($pre2, $indir, $ada2);

#get pairs
my $in1 = $pre1."_trimEndAda.fastq";
my $in2 = $pre2."_trimEndAda.fastq";

my $clean_pair1 = $pre1."_clean.fasta";
my $clean_pair2 = $pre2."_clean.fasta";
my $single_fasta  = $pre1."_".$pre2."_cleanSingle.fasta";
my $pair_stat = $pre1."_".$pre2."stat";

if ( !(-e $indir/$in1) or !(-e $indir/$in2) or (-e $indir/$clean_pair1) or 
		(-e $indir/$clean_pair2) or (-e $indir/$single_fasta) or  (-e $indir/$pair_stat) )
{die "!(-e $indir/$in1) or !(-e $indir/$in2) or (-e $indir/$clean_pair1) or (-e $indir/$clean_pair2) or (-e $indir/$single_fasta) or  (-e $indir/$pair_stat) "}

else 
{
	my $cmd_pair = "/Users/kaitang/Desktop/perl_code/Project_5/find_pairs_for_pair_end_reads_v2.pl $indir/$in1 $indir/$in2 $indir/$clean_pair1 $indir/$clean_pair2 $indir/$single_fasta > $indir/$pair_stat";
	
	print "$cmd_pair\n\n getting pair\n";
	`$cmd_pair`; 
}

#soap 

$soap_pe = "soap "







  exit;

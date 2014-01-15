#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 17 Nov 2010

like a .sh script to batch run BWA and samtools.
This is version 1, all use default parameter
=cut

my $debug = 0;

use strict;
my $usage = "$0 <input_dir> <out_dir>";
die $usage unless(@ARGV == 2);

my ($indir,$outdir) = @ARGV[0..1];

opendir (INDIR, $indir)
  or die "Cannot open dir $indir: $!";
  
 my @files = grep {/\.fast.$/} readdir INDIR;
 
 if ($debug)
 {
 	print STDERR join "\n",@files;	
 }
 
# if ($#files != 1)
#	{print STDERR join "\t",readdir INDIR;
#		die "not exact two files!";}

foreach my $file (@files)
{
	my $pre = "pre";
	if ($file =~ /(\S+)\.fast.$/)
	{$pre = $1;}
	my $bwa_out_sai = $pre."_n2_l25_k1.sai";
	
	my $bwa_aln_cmd = "time bwa aln -n 2 -k 1 -l 25 ~/Desktop/C24_genome/BWA_index/C24_TAIR9_5Chr.fas $indir/$file > $outdir/$bwa_out_sai";
	
	my $sam_file = $pre."_n2_l25_k1.sam";
	
	my $bwa_samse_cmd = "bwa samse  ~/Desktop/C24_genome/BWA_index/C24_TAIR9_5Chr.fas $outdir/$bwa_out_sai $indir/$file > $outdir/$sam_file";
	
 my $bam_file = "unname.bam";
 my $sam_pre = "sam_pre";
 
 if ($sam_file =~ /(\S+)\.sam/) {$sam_pre = $1;}
 $bam_file = $sam_pre.".bam";
 my $sam_to_bam_cmd = "samtools view -bS $outdir/$sam_file > $outdir/$bam_file";

my $sorted_bam = $sam_pre."_sorted";
 my $sort_cmd = "samtools sort $outdir/$bam_file $outdir/$sorted_bam";
 
 my $real_sorted_bam = $sam_pre."_sorted.bam";
 my $index_cmd = "samtools index $outdir/$real_sorted_bam";
 
 print "$bwa_aln_cmd\n\n\n\n";
 `$bwa_aln_cmd`;
 print "\n\n";
 
  print "$bwa_samse_cmd\n\n\n\n\n";
 `$bwa_samse_cmd`;
  print "\n\n";

  print "$sam_to_bam_cmd\n\n\n\n\n";
 `$sam_to_bam_cmd`;
  print "\n\n";

  print "$sort_cmd\n\n\n\n\n";
 `$sort_cmd`;
  print "\n\n";

  print "$index_cmd\n\n\n\n\n";
 `$index_cmd`;

}

exit;

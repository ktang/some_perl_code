#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 17 Nov 2010

like a .sh script to batch run BWA and samtools.
This is version 1, all use default parameter
=cut

use strict;
my $usage = "$0 <input_dir> <out_dir>";
die $usage unless(@ARGV == 2);

my ($indir,$outdir) = @ARGV[0..1];

opendir (INDIR, $indir)
  or die "Cannot open dir $indir: $!";
  
 my @files = grep {/\.fast.$/} readdir INDIR;
 
 if ($#files != 1)
	{print STDERR join "\t",readdir INDIR;
		die "not exact two files!";}

my ($pre1,$pre2) = ("pre1","pre2");	
	
if ($files[0] =~ /(\S+)\.fast.$/)
{$pre1 = $1;}

if ($files[1] =~ /(\S+)\.fast.$/)
{$pre2 = $1;}

my $bwa_out_1 = $pre1."_n2_k1_l25.sai";
	
 my $bwa_aln_cmd_1 = "time bwa aln -n 2 -k 1 -l 25 ~/Desktop/C24_genome/BWA_index/C24_TAIR9_5Chr.fas $indir/$files[0] > $outdir/$bwa_out_1";
 
 my $bwa_out_2 = $pre2."_n2_k1_l25.sai";
	
 my $bwa_aln_cmd_2 = "time bwa aln -n 2 -k 1 -l 25 ~/Desktop/C24_genome/BWA_index/C24_TAIR9_5Chr.fas $indir/$files[1] > $outdir/$bwa_out_2";
 
 my $sam_file = $pre1."_".$pre2."_n2_k1_l25_pe_default.sam";
 my $bwa_sampe_cmd = "bwa sampe ~/Desktop/C24_genome/BWA_index/C24_TAIR9_5Chr.fas $outdir/$bwa_out_1 $outdir/$bwa_out_2 $indir/$files[0] $indir/$files[1] > $outdir/$sam_file";
 
 my $bam_file = "unname.bam";
 my $sam_pre = "sam_pre";
 if ($sam_file =~ /(\S+)\.sam/) {$sam_pre = $1;}
 $bam_file = $sam_pre.".bam";
 my $sam_to_bam_cmd = "samtools view -bS $outdir/$sam_file > $outdir/$bam_file";
 
 my $sorted_bam = $sam_pre."_sorted";
 my $sort_cmd = "samtools sort $outdir/$bam_file $outdir/$sorted_bam";
 
 my $real_sorted_bam = $sam_pre."_sorted.bam";
 my $index_cmd = "samtools index $outdir/$real_sorted_bam";
 
 print "$bwa_aln_cmd_1\n\n\n\n";
 `$bwa_aln_cmd_1`;
 print "\n\n";
  print "$bwa_aln_cmd_2\n\n\n\n\n";
 `$bwa_aln_cmd_2`;
  print "\n\n";

  print "$bwa_sampe_cmd\n\n\n\n\n";
 `$bwa_sampe_cmd`;
  print "\n\n";

  print "$sam_to_bam_cmd\n\n\n\n\n";
 `$sam_to_bam_cmd`;
  print "\n\n";

  print "$sort_cmd\n\n\n\n\n";
 `$sort_cmd`;
  print "\n\n";

  print "$index_cmd\n\n\n\n\n";
 `$index_cmd`;

  exit;

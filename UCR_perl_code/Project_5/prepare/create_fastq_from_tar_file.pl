#! /usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1  Nov22, 2010

Nov 22,2010

this script is specified for tar file in
/Users/kaitang/Desktop/zhu_lab_data_orignal/100822_HWUSI-EAS1501_00007_61RF3AAXX_PE1

usage is
./create_fastq_from_tar_file.pl s_5_1_ /Users/kaitang/Desktop/zhu_lab_data_orignal/100822_HWUSI-EAS1501_00007_61RF3AAXX_PE1/Col_gl1/orignal_files /Users/kaitang/Desktop/zhu_lab_data_orignal/100822_HWUSI-EAS1501_00007_61RF3AAXX_PE1/Col_gl1/fastq/seq5/s_5_1_Col_gl1.fastq;
=cut
use strict;

my $debug = 0;
my $usage = "<$0> <pre_name> <indir> <outputfile>";

die $usage unless (@ARGV == 3);

my ($pre_name,$indir,$outfile) = @ARGV[0..2];

opendir (INDIR,"$indir")
	or die "cannot open indir $indir:$!";

open (OUT, ">$outfile")
	or die "cannot open output file $outfile:$!";

my @files = grep /$pre_name/, readdir INDIR;
if ($debug) {print STDERR join ("\t",@files),"\n"}
foreach my $file (sort @files)
{
	if ($debug) {print STDERR "$file\n";}
	open (IN ,"$indir/$file")
		or die "cannot open input_file $file:$!";
	my $line;
	while ($line = <IN>)
	{
		chomp $line;
		my @pts = split "\t",$line;
		if ($debug)
		{
			print STDERR join "\n",@pts;
		}
		my $p1 = "@".$pts[0];
		my $p2 = join ":",@pts[1..5];
		my $p3 = join "/",@pts[6..7];
		my @ls;
		$ls[0] = $p1."_".$p2."\#".$p3;
		$ls[1] = $pts[8];
		$ls[2] = "+";
		$ls[3] = $pts[9];
		print OUT join "\n",@ls;
		print OUT "\n";
	}
}

exit;

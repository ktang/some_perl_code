#!/usr/bin/perl -w
#original: /Users/tang58/deep_seq_analysis/Ros3_bind_RNA/figures/input
#this may be not updated.
use strict;

my $usage = "$0 <indir>";

die $usage unless (@ARGV == 1);

my $indir = $ARGV[0];

opendir(DIR, $indir) or die "cannot opendir $indir:$!";

my @files = grep /profile.txt$/, readdir(DIR);

foreach my $file(@files){
	if($file =~ /(\S+)_trim/){
		my $pre = $1;
		my $cmd_head = "head -2 $file > base_dist.txt";
		my $cmd_tail = "tail -21 $file > tail.txt";
		print $cmd_head, "\n";
		`$cmd_head`;
		print $cmd_tail, "\n";
		`$cmd_tail`;
		open (IN, $file) or die "cannot open $file:$!";
		open (OUT, ">length_dist.txt") or die "cannot open output:$!";
		<IN>; <IN>; <IN>;
		W : while (<IN>){
			if (/^$/){last W}
			print OUT $_;
		}
		
		
	}
}
#!/usr/bin/perl -w
#use hash to record the mapping situation

use strict;

my $usage = "$0 <dir> <file>";
die $usage unless(@ARGV == 2 );
my $indir = $ARGV[0];
my $infile = $ARGV[1];

#opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my $pre = "NONE";
if ($infile =~ /(\S+).soap$/){
	$pre = $1;
}

my $output = $pre."_profile_perl.txt";
if (-e "$indir/$output" or $pre eq "NONE") {
	die "$output exists or wrong name!!!\n";
}

open (IN,"$indir/$infile") or die "cannot open $infile";
open (OUT,">$indir/$output") or die "cannot opne $output";
#six counter.
my ($one,$two, $three, $more, $perfect_match)=
	(0,0,0,0,0);
	
my (%trans, %genes);
my %perfect_flag;
my %marks;

my $last = "NONE";
my $gene = "NO";
my $no_mismatch_flag = 0;

my $i = 0;

while(<IN>){
	chomp;
	$i++;
	my @a = split "\t";
	if ($a[7] =~ /(AT.G.....)\.(\d+)/){
		my $gene = $1;
		my $version = $2;
		$marks{$a[0]}->{$gene}->{$version}++;
		if ($a[9] == 0){
			$perfect_flag{$a[0]} = 1;
		}
	}
	else{
		print STDERR "wrong gene id:$_\n";
	}
}
close(IN);

foreach my $read (keys %marks){
	my $hits = scalar (keys %{$marks{$read}});
	
	#my $max = 0;
	foreach my $gene (keys %{$marks{$read}}){
		my $max = 1;
		foreach my $ver (keys %{$marks{$read}->{$gene} }){
			if ($marks{$read}->{$gene}->{$ver} > $max ){
				$max = $marks{$read}->{$gene}->{$ver};
			}
		}
		if ($max > 1){
			$hits += ($max - 1);
		}
	}
	
	if ($hits == 1){
			$one ++;
			if (defined $perfect_flag{$read}){
				$perfect_match++;
			}
		}elsif($hits == 2){
			$two++;
		}elsif($hits == 3){
			$three++;
		}elsif($hits > 3){
			$more++;
		}else{
			
				print STDERR "wrong hit number: $hits\n$_\n";
			
		}
	}
print OUT join("\t",("library",">3_location","3_loacation","2_loc","1_loc" ,"perfect_match", )), "\n";
print OUT join("\t", ($pre,$more, $three, $two, $one, $perfect_match)), "\n";
close(OUT);

for my $key (keys %trans){
	print STDERR $key,"\n";
}

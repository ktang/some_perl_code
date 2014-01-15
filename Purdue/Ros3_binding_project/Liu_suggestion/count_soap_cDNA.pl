#!/usr/bin/perl -w
#this version may be wrong
#as do not consider the situation 
#that one read can have two or more
#hits on the same transcript.
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
	(0,0,0,0,0,0);
	
my (%trans, %genes);

my $last = "NONE";
my $gene = "NO";
my $no_mismatch_flag = 0;

my $i = 0;

while(<IN>){
	chomp;
	$i++;
	my @a = split "\t";
	if ($a[0] ne $last){#new start
		my $hits = scalar (keys %genes);
		%genes = ();
		undef %genes;
		$gene = "NO";
		%trans = ();
		undef %trans;
		
		if (defined $trans{$a[7]}){
			print STDERR "same transcript have two hits:\n$_\n";
		}else{
			$trans{$a[7]} = 1;
		}
		
		if ($a[7] =~ /(AT.G.....)\./){
			$gene = $1;
		}else{
			die "wrong gene id:$_\n"
		}
		
		if ($gene ne "NO"){
			$genes{$gene} = 1;
		}else{
			die "wrong gene $gene\n"
		}
		
		if ($hits == 1){
			$one ++;
			if ($no_mismatch_flag == 1){
				$perfect_match++;
			}
		}elsif($hits == 2){
			$two++;
		}elsif($hits == 3){
			$three++;
		}elsif($hits > 3){
			$more++;
		}else{
			if($i != 1){
				print STDERR "wrong hit number: $hits\n$_\n";
			}	
		}
		if ($a[9] == 0){
			$no_mismatch_flag = 1;
		}else{
			$no_mismatch_flag = 0;
		}
		
		$last = $a[0];
	}else{
		if (defined $trans{$a[7]}){
			print STDERR "same transcript have two hits:\n$_\n";
		}else{
			$trans{$a[7]} = 1;
		}
		
		if ($a[7] =~ /(AT.G.....)\./){
			$gene = $1;
		}else{
			die "wrong gene id:$_\n"
		}
		
		if ($gene ne "NO"){
			$genes{$gene} = 1;
		}else{
			die "wrong gene $gene\n"
		}
		
		if ($a[9] == 0){$no_mismatch_flag = 1}
		
	}
}

close(IN);
my $hits = scalar (keys %genes);
if ($hits == 1){
	$one ++;
	if ($no_mismatch_flag == 1){
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

print OUT join("\t",("library",">3_location","3_loacation","2_loc","1_loc" ,"perfect_match", )), "\n";
print OUT join("\t", ($pre,$more, $three, $two, $one, $perfect_match)), "\n";
close(OUT);

for my $key (keys %trans){
	print STDERR $key,"\n";
}

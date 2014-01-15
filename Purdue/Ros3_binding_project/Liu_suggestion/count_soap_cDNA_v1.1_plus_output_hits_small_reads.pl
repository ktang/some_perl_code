#!/usr/bin/perl -w
#use hash to record the mapping situation

use strict;

my $debug = 0;

my $usage = "$0 <dir> <file>";
die $usage unless(@ARGV == 2 );
my $indir = $ARGV[0];
my $infile = $ARGV[1];

#opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my $pre = "NONE";
if ($infile =~ /(\S+).soapChr$/){
	$pre = $1;
}

my $output = $pre."_cDNA_profile.txt";
my $soapout = $pre."_cDNA_back2Chr.soap2";
if (-e "$indir/$output" or $pre eq "NONE" or -e "$indir/$soapout") {
	die "$output exists or wrong name!!!\n";
}

open (IN,"$indir/$infile") or die "cannot open $infile";
open (OUT,">$indir/$output") or die "cannot open $output";
open (OUT2,">$indir/$soapout") or die "cannot open $soapout";

my ($one,$two, $three, $more, $perfect_match)=
	(0,0,0,0,0);
	
my %genes;
my %perfect_flag;
my %marks;

#my $gene = "NO";
my $no_mismatch_flag = 0;

#my $i = 0;

my @ins = <IN>;
close(IN);

#while(<IN>){
for (my $i = 0; $i <= $#ins; $i++){
	my $this = $ins[$i];
	chomp $this;
	my @a = split "\t", $this;
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


my %outs; #record hits number

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
#		if ($max > 1){
			$hits += ($max - 1);
#		}
	}
	
	if ($hits == 1){
			$one ++;
			$outs{$read} = 1;
			if (defined $perfect_flag{$read}){
				$perfect_match++;
			}
	}
	
	elsif($hits == 2){
			$outs{$read} = 2;
			$two++;
	}
	
	elsif($hits == 3){
			$outs{$read} = 3;
			$three++;
	}
	
	elsif($hits > 3){
			$more++;
	}
	
	else{
			print STDERR "wrong hit number: $hits\n$_\n";
	}
}

print OUT join("\t",("library",">3_location","3_loacation","2_loc","1_loc" ,"perfect_match", )), "\n";
print OUT join("\t", ($pre,$more, $three, $two, $one, $perfect_match)), "\n";
close(OUT);

my $last = "NONE";

my (%gene_out, %trans_out);

#for my $key (keys %outs){
#	print STDERR $outs{$key}, "\n";
#}

for (my $i = 0; $i <= $#ins; $i++){
	my $this = $ins[$i];
	chomp $this;
	my @a = split "\t", $this;
	
	if (defined $outs{$a[0]}){
		
	#		my $hits = $outs{$a[0]}; 
	#		print STDERR $hits,"\n";
		
		if($a[0] ne $last){
			$last = $a[0];
			
			if( $outs{$a[0]} == 1){
				print OUT2 out_line($this, 1), "\n";
			}
			
			elsif  ($outs{$a[0]} == 2){
				$trans_out{$a[0]} -> {$a[7]} = 1;
				if($a[7] =~ /(AT.G.....)\./) {
					$gene_out{$a[0]} ->{$1} = 1;
					if ($debug) { print STDERR $1, "	:if_2\n"}
				}
				print OUT2 out_line($this, 2), "\n";
			}
			
			elsif ($outs{$a[0]} == 3){
				print OUT2 out_line($this, 3), "\n";
				$trans_out{$a[0]} -> {$a[7]} = 1;
				if($a[7] =~ /(AT.G.....)\./) {
					$gene_out{$a[0]} ->{$1} = 1;
					if ($debug) { print STDERR $1, "	:if_3\n"}
				}
			}
		}
		
		else{
			
			if($outs{$a[0]}  == 1){
				next;
			}
			
			if ($outs{$a[0]}  == 2){
				my $this_transcript = $a[7];
				my $this_gene;
				if($a[7] =~ /(AT.G.....)\./) {
					$this_gene = $1;
					if ($debug) { print STDERR $1, "	:else_2\n"}
				}
				
				die "this_gene not defined" unless (defined $this_gene);
				
				if (defined $trans_out{$a[0]}->{$this_transcript} or !(defined $gene_out{$a[0]}->{$this_gene}) ){
					print OUT2 out_line($this, 2), "\n";
					$gene_out{$a[0]}->{$this_gene} = 1;
				}
			}
			
			if ($outs{$a[0]}  == 3){
				my $this_transcript = $a[7];
				my $this_gene;
				if($a[7] =~ /(AT.G.....)/) {
					$this_gene = $1;
					if ($debug) { print STDERR $1, "	:else_3\n"}
				}
				
				die "this_gene not defined" unless (defined $this_gene);
				
				if (defined $trans_out{$a[0]}->{$this_transcript} or !(defined $gene_out{$a[0]}->{$this_gene})  ){
					print OUT2 out_line($this, 3), "\n";
					$gene_out{$a[0]}->{$this_gene} = 1;
				}			
			}
		}
	}
}

sub out_line{
	my ($line, $hit) = ($_[0], $_[1]);
	
	my @a = split "\t", $line;
	$a[3] = $hit;
	my $gene_cDNA = $a[7];
	my $start_cDNA = $a[8];
	my $strand_cDNA = $a[6];
	my ($chr, $start, $optional) = split "_", $a[-1];
	if (defined $optional){
		$start = $optional;
		if($a[6] eq "+"){$a[6] = "-"}
		elsif($a[6] eq "-") {$a[6] = "+"}
		else {die "wrong strand symble!!"}
	}
	$a[10] = join("_",($gene_cDNA, $start_cDNA ,$strand_cDNA));
	$a[7] = $chr;
	$a[8] = $start;
	return join("\t", @a[0..10]);
}

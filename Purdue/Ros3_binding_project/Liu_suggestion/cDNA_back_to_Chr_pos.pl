#!/usr/bin/perl -w

use strict;

my $debug = 0;

my $usage = "$0 <dir> <file>";
die $usage unless(@ARGV == 2 );
my $indir = $ARGV[0];
my $infile = $ARGV[1];



my $gff_file = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_GFF3_genes_transposons.gff";
if ($debug){
	$gff_file = "/Users/tang58/try/Ros3_binding/batch_soap/cDNA_back_to_genome/gff.txt";
}

open (GFF, $gff_file) or die "cannot open $gff_file";
#my $gffs = <GFF>;
#close (GFF);

my (%strands, %exons);

while (<GFF>){
	chomp;
	my @a = split "\t";
	if ($a[2] =~ /exon/){
		my ($parent,$id) = split "=",$a[8];
		$strands{$id} = $a[6];
#		my $start = $a[3];
#		my $length = $a[4] - $a[3] + 1;
		if ($a[6] eq "+"){
			my $start =  $a[3];
			my $length = $a[4] - $a[3] + 1;
			push @{$exons{$id}}, [$start, $length];
		}elsif ($a[6] eq "-"){
			my $end = $a[4];
			my $length = $a[4] - $a[3] + 1;
			push @{$exons{$id}}, [$end, $length];
		}else{
			die "wrong strand symble:$a[6]\n\n";
		}
	}else{
		next;
	}
}

#opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my $pre = "NONE";
if ($infile =~ /(\S+).soap$/){
	$pre = $1;
}

my $output = $pre."_chr_info.soap";
if (-e "$indir/$output" or $pre eq "NONE") {
	die "$output exists or wrong name!!!\n";
}

open (IN,"$indir/$infile") or die "cannot open $infile";
open (OUT,">$indir/$output") or die "cannot opne $output";

while (<IN>){
	chomp;
	my @a = split "\t";
	my $id = $a[7];
	if(defined $strands{$id} and defined $exons{$id}){
		my $pos_trans = $a[8];
		my @exon_array = @{$exons{$id}};
		my $left = $a[8];
		
		if ($strands{$id} eq "+"){
			for (my $i = 0; $i <= $#exon_array; $i++){
				my ($start, $length) = @{$exon_array[$i]};
				if ($left <= $length){
					
					my $chr_pos = $start + $left -1;
					my $Chr = "NONE";
					if ($id =~ /AT(.)G/){
						$Chr = "Chr".$1;
					}else{die "wrong chr:	$id"}
					my $last_col = join ("_", ($Chr, $chr_pos));
					print OUT join ("\t", (@a, $last_col)), "\n";
					last;
				}
				
				else{
					$left -= $length;
					if ($left < 0){
						print STDERR "exceed the gene length\n";
					}
				}
			}
		}
		
		elsif($strands{$id} eq "-"){
			for (my $i = 0; $i <= $#exon_array; $i++){
				my ($end, $length) = @{$exon_array[$i]};
				
				if ($left <= $length){
					
					my $chr_pos = $end - $left + 1;
					my $Chr = "NONE";
					if ($id =~ /AT(.)G/){
						$Chr = "Chr".$1;
					}else{die "wrong chr:	$id"}
					
					my $chr_start = $chr_pos - $a[5] + 1;
					my $last_col = join ("_", ($Chr, $chr_pos, $chr_start));
					print OUT join ("\t", (@a, $last_col)), "\n";
					last;
				}
				else{
					$left -= $length;
					if ($left < 0){
						print STDERR "exceed the gene length\n";
					}
				}
			}
		}
		
		else{
			die "wrong strand in strands{$id}:$strands{$id}";
		}
	}
	
	else{
		print STDERR "no $id in gff\n";
	}
}#while
close(IN);
close(OUT);

if ($debug){
	foreach my $id (sort keys %exons){
		print STDERR "$id:";
		my @array = @{$exons{$id}};
		for (my $i = 0; $i <= $#array; $i++){
				my ($start, $length) = @{$array[$i]};	
				print STDERR "$start,$length\n";
		}
	}
}
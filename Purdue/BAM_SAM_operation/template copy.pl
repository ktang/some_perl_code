#!/usr/bin/perl -w
use strict;
use File::Spec;

#v0.0
# start or end is in structural RNA

#v0.1
#whole read is in structural RNA

my $debug = 0;

my $gff_file = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/TAIR10_GFF3_structural_RNA_simpleChr.txt";
die unless (-e $gff_file);

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <outdir> <outpre>\n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;

die unless (-e $input) ;
die unless (-d $outdir);

my $out_no_structural_RNA	= File::Spec->catfile ($outdir, $outpre . "_without_structural_RNA.bam");
my $out_with_structural_RNA	= File::Spec->catfile ($outdir, $outpre . "_with_structural_RNA.bam");
my $out_nomap			= File::Spec->catfile ($outdir, $outpre . "_nomap.bam");
my $out_stat			= File::Spec->catfile ($outdir, $outpre . "_length_stat.txt");

die if(-e $out_no_structural_RNA) ;
die if(-e $out_with_structural_RNA) ;
die if(-e $out_nomap);
die if(-e $out_stat);

my @category = ("no_structural", "structural", "no_hit");

my %length_dist;

if ($debug){
	print STDOUT join("\n", ( $out_no_structural_RNA, $out_with_structural_RNA, $out_nomap, $out_stat )), "\n\n";
	exit;
}

my %h;
record ( $gff_file , \%h);

open( IN, "samtools view -h $input |") or die;
open( NS ,  "| samtools view -b -S -  > $out_no_structural_RNA") or die;
open( WS ,  "| samtools view -b -S -  > $out_with_structural_RNA") or die;
open( NO ,  "| samtools view -b -S -  > $out_nomap") or die;

open(STAT, ">>$out_stat") or die;

my $last_id = "NA";
my $num_of_hit = 0;
my @table = ();
my $flag_no_SR = 1 ;

# structural_RNA
while(<IN>){
	if(/^@/){
		print NS $_;
		print WS $_;
		print NO $_;
		#print  $_;
	}else{
		chomp;
		my @temp = split "\t";
		my ($chr, $pos ) = ($temp[2], $temp[3]);
		my $length = length ( $temp[9]);
		my $end = $pos + $length - 1;
		if ($chr eq "*"){
			print NO $_, "\n";
			my $len = length ( $temp[9]);
			$length_dist{"no_hit"}->{$len}++;
			next;
		}
		if( $temp[0] ne $last_id ){
			
			if( $last_id ne "NA" ){
				if($flag_no_SR){
					my $len = length ( $table[0]->[9] );
					$length_dist{"no_structural"}->{$len}++;
					for my $i(0..$#table){
						print NS join("\t", @{$table[$i]}, "NH:i:" . $num_of_hit), "\n"; #NH:i:1
					}
				}else{
					my $len = length ( $table[0]->[9] );
					$length_dist{"structural"}->{$len}++;
					for my $i(0..$#table){
						print WS join("\t", @{$table[$i]}, "NH:i:" . $num_of_hit), "\n"; #NH:i:1
					}
				}
			}
			
			$last_id = $temp[0];
			$flag_no_SR = 1;
			$num_of_hit = 1;
			@table = ();
			push @table, [@temp];
			
			if ( defined  $h{$chr}->[$pos] and  defined  $h{$chr}->[$end] ){
				$flag_no_SR = 0;
			}
		}else{
			$num_of_hit++;
			push @table, [@temp];
			if ( defined  $h{$chr}->[$pos]  and  defined  $h{$chr}->[$end] ){
				$flag_no_SR = 0;
			}
		}
	}
}

#########
#	last one
############
if($flag_no_SR){
	my $len = length ( $table[0]->[9] );
	$length_dist{"no_structural"}->{$len}++;
	for my $i(0..$#table){
		print NS join("\t", @{$table[$i]}, "NH:i:" . $num_of_hit), "\n"; #NH:i:1
	}
}else{
	my $len = length ( $table[0]->[9] );
	$length_dist{"structural"}->{$len}++;
	for my $i(0..$#table){
		print WS join("\t", @{$table[$i]}, "NH:i:" . $num_of_hit), "\n"; #NH:i:1
	}
}

close IN;
close WS;
close NS;
close NO;

for my $class (@category){
	print  STAT  $class, "\n";
	my $total = 0;
	foreach my $len (sort {$a <=> $b} keys %{$length_dist{$class}}){
		print STAT  join("\t", ( $len,  $length_dist{$class}->{$len})) , "\n";
		$total+= $length_dist{$class}->{$len};
	}
	print STAT  join("\t", ( "total",  $total)) , "\n";
	print  STAT "\n";
}

close STAT ;
exit;

#0	  1	 2	 3	 4	 5	 6	 7	 8
#1       TAIR10  gene    111890  111961  .       -       .       ID=AT1G01270;Note=tRNA;Name=AT1G01270
#1       TAIR10  gene    306384  306456  .       +       .       ID=AT1G01870;Note=tRNA;Name=AT1G01870
#record ( $gff_file , \%h);
sub record{
	my ($file, $ref) = @_;
	die unless (-e $file);
	open(SUB, $file) or die;
	
	while(<SUB>){
		chomp;
		my @a = split "\t";
		for my $i($a[3] .. $a[4]){
			$ref->{$a[0]}->[$i] = 1;
		}
	}
	
	close SUB;	
}
#!/usr/bin/perl -w
use strict;
use File::Spec;

# v0.2
# modified from v0.1
# and use samtools depth not mpileup
# debug in /Users/tang58/misc/Huiming_Zhang/RDM15_At4g31880/sRNA_20130929/DSR


#v0.1
#the raw srcipt use rrp6l1 as mut
#in the v0.1 version, the mutant bam file is as input


#v0.0
#step 4 reduce boundary to the first and last bp of the two combind file

my $debug = 0;

my $run_debug = 0;



my %format_h = ( "bed" => 1,
		"coordinate" => 1);


# 8/2 = 4
# 8:dividend
#2: divisor
#my $usage = "$0 \n <input_list> <format (bed or coordinate)> <outdir> <outpre> \n\n";
#die $usage unless(@ARGV == 4);

my $usage = "$0 \n <input_list> <format (bed or coordinate)> <outdir> <outpre> <bam_file1> [ <bam2> ...] \n\n";
die $usage unless(@ARGV >= 5);

my $input = shift or die;
my $format = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;

die unless (@ARGV>=1);

my @bams;
for my $tmp(0..$#ARGV){
	my $f = shift or die;
	die unless ( -e $f);
	push @bams, $f;
}



die unless ( -e $input);
#die if( -e $output );
die "coordinate  bed" unless ( defined $format_h{$format});
#die unless (defined $labels{$label});

#my $output = File::Spec->catfile($outdir, $outpre . "_reduced_boundary.txt");
my $output = File::Spec->catfile($outdir, $outpre . "_RB.txt");
die "$output exists!!\n\n" if( -e $output );



if($debug){
	print STDERR "bams:\n";
	print STDERR join("\n", @bams), "\n\n";
	exit;
}

#open( IN, "samtools view $input |") or die;
#open( , ">>$output") or die "cannot open $output: $!";

my @list;

if($format eq "bed"){
	read_bed( $input, \@list );
}elsif( $format eq "coordinate" ){
	read_coor( $input, \@list );
}
else{
	die $format;
}

open(OUT, ">>$output") or die;

my @head_h = split "\t", $list[0];

print OUT join("\t", ("coor_reduced", "coor_raw", @head_h[1..$#head_h] ) ), "\n";
for my $i (1..$#list){
	my $this = $list[$i];
	my @a = split "\t", $this;
	my $coor_raw = $a[0];
	
	if ( $run_debug){
		print STDERR $this, "\n";
	}
	
	my ($chr, $s, $e) = get_coor($coor_raw);
	
	if($s == -1){
		print $this, "\n\n";
		die;
	}
	
#	my $first_line = `samtools mpileup -f $fa_h{$chr} -r $coor_raw $WT_bam $rrp6_bam 2> /dev/null | head -1 | awk '{print \$2}' `;
#	my $last_line  = `samtools mpileup -f $fa_h{$chr} -r $coor_raw $WT_bam $rrp6_bam 2> /dev/null | tail -1 | awk '{print \$2}' `;
	my ($first_line, $last_line) = (-1, -1);
	my @tmp = ();
	@tmp = `samtools depth -r $coor_raw  @bams | awk '{print \$2}' `;
	$first_line = $tmp[0];
	$last_line  = $tmp[-1];
	chomp( $first_line );
	chomp( $last_line );
	
	print STDERR $coor_raw, "\n" if ($first_line >= $last_line);
	
	if($run_debug){
		print STDERR "first:", $first_line, "\n";
		print STDERR "last:", $last_line,  "\n";
	}
	
	my $reduced_coor = "NA";
	if($first_line =~ /\d+/){
		$reduced_coor = "$chr:$first_line-$last_line";
		print OUT join("\t",$reduced_coor , @a), "\n";
	}else{
		print $a[0], "\n";
	}
	
	
	
}
close OUT;
exit;

sub read_bed{#( $input, \@list );
	my ($file, $ref) = @_;
	die unless( -e $file);
	open(IN, $file) or die;
	
	while(<IN>){
		next if(/^browser/);
		next if(/^track/);
		next if(/^#/);
		next unless(/\w+/);
		last;
	}
	
#	my $h = <IN>;
	my $h = $_;
	chomp $h;
	my @a_h = split "\t", $h;
	
	my $last_index = $#a_h;
	
	$ref->[0] = join("\t", ("coordinate", @a_h[3..$last_index]));
	my $i = 0;
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		die if ($a[0] =~ /:/);
		my $chr = simple_chr ($a[0] );
		my $coor = "$chr:$a[1]-$a[2]";
		$ref->[$i] = join("\t", ( $coor, @a[3..$last_index]));
	}
	close IN;
}
#( $input, \@list  );
sub read_coor{ 
	my ($file, $ref) = @_;
	die unless( -e $file);
	open(IN, $file) or die;
	
	while(<IN>){
		next if(/^browser/);
		next if(/^track/);
		next if(/^#/);
		next unless(/\w+/);
		last;
	}
	
	if($run_debug){
		print STDERR $file, "\n";
	}
#	my $h = <IN>;
	my $h = $_;
	chomp $h;
	$ref->[0] = $h;
	my $i = 0;
	
	while(<IN>){
		$i++;
		chomp;
		my @a = split "\t";
		die  $a[0] unless ($a[0] =~ /:/);
	#	my $chr = simple_chr ($a[0] );
	#	my $coor = "$chr:$a[1]-$a[2]";
		if( $run_debug ){
			print STDERR "read_coor: ", join("\n",  @a ), "\n";
		}
		$ref->[$i] = join("\t", @a );
	}
	close IN;
}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
 #get_coor($coor_raw);
 sub get_coor{
	my ($coor) = @_;
	if($coor =~ /(\S+):(\d+)-(\d+)/){
		return ($1, $2, $3);
	}else{
		print STDERR $coor, "\n\n";
		return( -1 , -1 , -1);
	}
 }

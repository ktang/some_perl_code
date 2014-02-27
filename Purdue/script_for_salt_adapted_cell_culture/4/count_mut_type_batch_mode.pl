#!/usr/bin/perl -w

#/Users/tang58/scripts_all/perl_code/Purdue/script_for_salt_adapted_cell_culture/4/count_mut_type_batch_mode.pl out_anno/ > dep8_per30_anno_type_num.txt


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;


use strict;
use File::Spec;
my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> STDOUT\n\n";
die $usage unless(@ARGV == 1);


my $indir = shift or die "indir";
#my $outdir = shift or die "outdir";

die unless (-d $indir);
#die unless (-d $outdir);

opendir (DIR, $indir) or die;
my @files = grep /\.txt$/ , readdir DIR;
closedir DIR;

print STDERR "files: \n";
print STDERR join("\n", @files ), "\n\n";

#my @anno_types = qw/CDS intron UTR TE ncRNA pseudogene IG/;
my @anno_types = qw/CDS intron UTR TE ncRNA pseudoG IG/;
#  1 => "CDS",
#  2 => "TE",
#  3 => "ncRNA",
#  4 => "UTR",
#  5 => "intron",
#  6 => "pseudogene"
#  7 => "IG"
=head
foreach my $file (@files){
	my $pre;
	#baseAnno.txt
	#if ($file =~ /(\S+)_SelfPer/){
	if ($file =~ /(\S+)_baseAnno/){
		$pre = $1;
		push @labels, $pre;
	}else{
		die $file;
	}
}
=cut
my %nums; #h{mut_type}->file = 
#chr	pos	ref	mut_nrpe1_150	dep_nrpe1_150	per_nrpe1_150	seq_nrpe1_150	qual_nrpe1_150
#chr1	29081	G	G=>DEL	14	14.29	...,,,.,,,,-2ta,-2ta.,	DFFI)@J>##D4H<

#my $col_num = 4;

if ($debug) {
	exit;#code
}

my @labels= ();

foreach my $file (@files){
	my $pre;
	
	#if ($file =~ /(\S+)\.txt$/){
	if ($file =~ /(\S+)_baseAnno\.txt$/){
		$pre = $1;
		push @labels, $pre;
	}else{
		die $file;
	}
	
	my $input = File::Spec->catfile($indir, $file);
	
	open(IN, $input) or die "$input: $!";
	my $head = <IN>;
	while (<IN>) {
		chomp;
		my @a = split "\t";
		
		my $anno = $a[-2];
		if ( $anno eq "pseudogene") {
			$anno = "pseudoG" ;
		}
		
		
		my $this = $a[-3];
		my @paires = split ";", $this;
		if ( @paires==4 ) {
			if( $this =~ /DEP=(\d+);TYPE=(\S+);ALT=(\S+);PER=(\S+)$/ ){
				my ($dep, $type, $alt , $per ) = ($1, $2, $3, $4);
				
				my @types = split ",", $type;
				my @alts = split ",", $alt;
				my @pers = split ",", $per;
				my $ind_of_maximum = Kai_Module::get_ind_of_maximum(\@pers);
				my $max_per = $pers[$ind_of_maximum];
				
				my $type_max = $types[$ind_of_maximum];
				
				$nums{$type_max}->{$pre}->{$anno}++;
				
				
			}else{
				die $_;
			}
		}
		else{
			die $_;
		}
		
		#my $type = $a[ $col_num  - 1];
		#my $type = $a[ -2 ];
		#$nums{$type}->{$file} += 1;
	}
	
	close IN;
}

my @mut_types = qw/SNP INS DEL/;

foreach my $mut_type(@mut_types){
	print $mut_type, "\n";
	print join("\t", ("pos_type", @labels)), "\n";
	foreach my $anno_type (@anno_types){
		print $anno_type, "\t";
		foreach my $i(0..$#labels){
			my $n = 0;
			my $lab = $labels[$i];
			if (defined $nums{$mut_type}->{$lab}->{$anno_type} ) {
				$n = $nums{$mut_type}->{$lab}->{$anno_type};#code
			}
			print $n;
			if ($i == $#labels) {
				print "\n";
			}else{
				print "\t";
			}
			
			
		}
	}
	print "\n\n";
}

exit;
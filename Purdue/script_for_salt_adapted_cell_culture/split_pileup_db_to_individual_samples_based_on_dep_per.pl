#!/usr/bin/perl -w

#split_pileup_db_to_individual_samples_based_on_dep_per.pl
# input file is like
# Jan8_var_db_pos_14240loci_in_19bam_5MutBam_pileup_NoN_13826loci_left_categorize_NoMut_in5mut_1967loci.txt
# and
# pileupdb Jan8_var_db_pos_14240loci_in_19bam_pileup.txt
# use depth and per as cutoff to get loci information for each sample
#

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use strict;
use File::Spec;

my $debug = 1;

my $depth_cutoff = 8;
my $per_cutoff = 30;


#my $usage = "\n $0 \n <input_pileup_db> <input_categorize_file> <outdir> <outpre> <sample_number> \n\n";
#die $usage unless (@ARGV == 5);

my $usage = "\n $0 \n <input_pileup_db> <input_categorize_file> <outdir> <sample_number> \n\n";
die $usage unless (@ARGV == 4);

my $input_pileup_db = shift or die;
my $input_categorize_file = shift or die;
my $outdir = shift or die;

#my $outpre = shift or die;
my $sample_number = shift or die;

#my $postfix = shift or die;
#my $cutoff = shift or die;

#die unless (-e $input);
die unless (-e $input_pileup_db);
die unless (-e $input_categorize_file);
#die if (-e $output);
#open (OUT, ">$output") or die;
#open(IN, $input) or die;

my %pileup_db_h ; # pileup_db_h{head/join("_", @[0..2])}->{1..n} = join("\t", 3col)
read_pileup_db($input_pileup_db, \%pileup_db_h, $sample_number );

#make_output_files( $input_categorize_file, $outdir,  );
open(IN, $input_categorize_file) or die;
my $h = <IN>;
chomp $h;
my @h_a = split "\t", $h;
#print join("\t", (@h_a[0..2], "percentage", @h_a[3..$#h_a]) );

#my @fhr;
#my @fhw;
#for my $i(0..$last_index){
#	open($fhr[$i], "<" , $files[$i]) or die "$files[$i]";
#	open($fhw[$i], ">" , $outputs[$i]) or die $outputs[$i], ": $!";
#	print {$fhw[$i]} $head, "\n";
#}


my @output_file_names;
my @fhw;
#0	1
#33	34

for my $i (1..$sample_number){
	my $this = $h_a[ $i + 2 ];
	my $output_tmp = File::Spec->catfile( $outdir, $this . "_dep" .  $depth_cutoff . "_per" . $per_cutoff . ".txt"  );
	die if (-e $output_tmp);
	open( $fhw[$i], ">>$output_tmp" ) or die;
	print {$fhw[$i]} join("\t", (@h_a[0..2], $pileup_db_h{head}->{$i}, $this)) ,"\n";
}

while (<IN>) {
	chomp;
	my @a = split "\t";
	for my $i ( 1..$sample_number){
		my $this = $a[$i + 2];
		my @paires = split ";", $this;
		if ( @paires==4 ) {
			if( $this =~ /DEP=(\d+);TYPE=(\S+);ALT=(\S+);PER=(\S+)$/ ){
				my ($dep, $type, $alt , $per ) = ($1, $2, $3, $4);
				
				my @types = split ",", $type;
				my @alts = split ",", $alt;
				my @pers = split ",", $per;
				my $ind_of_maximum = Kai_Module::get_ind_of_maximum(\@pers);
				my $max_per = $pers[$ind_of_maximum];
				
				if ($dep >= $depth_cutoff and $max_per >= $per_cutoff) {
					my $tmp_lab = join("_", @a[0..2]);
					print  {$fhw[$i]} join("\t", (@a[0..2], $pileup_db_h{$tmp_lab}->{$i}, $this)) ,"\n";
				}
				
				
			}else{
				die $_;
			}
		}
	}
}



close IN;
exit;

#read_pileup_db($input_pileup_db, \%pileup_db_h);
#my %pileup_db_h ; # pileup_db_h{head/join("_", @[0..2])}->{1..n} = join("\t", 3col)
sub read_pileup_db{
	my ($file, $ref, $sample_number_sub) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	chomp $head;
	my @h_a_sub = split "\t", $head;
	for my $i ( 1..$sample_number_sub){
		$ref->{head}->{$i} = join("\t", ( @h_a_sub[($i * 3)..($i * 3 + 2)] ));
	}
	while ( <IN> ) {
		chomp;
		my @a = split "\t";
		my $tmp_lab = join("_", @a[0..2]); 
		for my $i ( 1..$sample_number_sub){
			$ref->{$tmp_lab}->{$i} = join("\t", ( @a[($i * 3)..($i * 3 + 2)] ));
		}
	}
	close IN;
}

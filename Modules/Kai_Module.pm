package Kai_Module;


use strict;
use warnings;
use File::Spec;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.0.0;
@ISA         = qw(Exporter);
@EXPORT      = qw(simple_chr round);
@EXPORT_OK   = qw(simple_chr round);
%EXPORT_TAGS = ( DEFAULT => [qw(&simple_chr)],
                 Both    => [qw(&simple_chr &round)]);
#%EXPORT_TAGS = ( DEFAULT => [qw(&func1)],
 #                Both    => [qw(&func1 &func2)]);



#my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);
#		$res_ref->[$i] = eval sprintf ("%.2f", 100 * ($numerator_ref->[$i] ) / ( $denominator_ref->[$i] ) );


=head
sub first_call{
	my ($x ) = @_;
	second_call($x);
}
sub second_call{
	my ($x) = @_;
	print STDERR $x, "\n\n";
}
=cut

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


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub get_sum{
	my @tmp = @_;
	my $s = 0;
	foreach my $i(@tmp){
		$s+=$i;
	}
	return $s;
}

sub add_chr{
	my ($chr) = @_;
	if ($chr =~ /^\d+$/) {
		return "Chr" . $chr;
	}elsif($chr eq "Mt" ){
		return "mitochondria";
	}elsif( $chr eq "Pt"){
		return "chloroplast"
	}else{
		die "wrong chr from add_chr: $chr\n\n";
	}
	
	#chloroplast, mitochondria
}

# split_pileup( $dep, $pileup_str, $qual_str )
#+	43
#-	45
##	=35
#$	36
#^	94
# as I require MAPQ>= 20, so there cannnot be ^+ and ^-
# maybe ^^ (seems max MAPQ = 42)

#  56789:;<=>?@ABCDEFGHIJK

=head
samtools mpileup  -f /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome.fa
-r 1:973010-973023 N11_S20_XZ7_ddc_150mM_bowtie2_XH_Dec24_sort_rmdup.bam 2>/dev/null |  \less

1       973010  C       15      ,$,,....,..,..,.        EH?EGGG*JI.JF>F
1       973011  A       14      ,,....,-2at..,..,.      HCGIEJ0JI=JG=F
1       973012  A       14      ,,....*..c..,.  H?G3HF!J9)JI=F
1       973013  T       14      ,,...-2TA.-4TATA*.-2TA.-2TA,-2ta.-2TA.-2TA,-2ta.-2TA    HCI=DA!DD?DD@D
1       973014  T       14      ,,..**,*******  HDI=!!!!!!!!!!
1       973015  A       14      ,,..**t*******  F>FC!!!!!!!!!!
1       973016  T       14      ,,...*,..,..,.  DDDD!!!!!!!!!!
1       973017  A       14      ,,..T*,..,..,.  CCCC!!!!!!!!!!
1       973018  T       14      ,,....,..,..,.  BBB@!!!!!!!!!!

, , . . .-2TA	.-4TATA	* .-2TA	.-2TA ,-2ta	.-2TA	.-2TA	,-2ta	.-2TA
H C I = D	A	! D	D	?	D	D	@	D
=cut

sub split_pileup{
	my $debug = 0;
	my $print_results = 0;
	
	my ( $dep, $pileup_str, $qual_str ) = @_;
	
	if ($dep == 0) {
		return "";
	}
	
	
	my @quals = split "", $qual_str;
		
	my @pileups = ();
	my $pileup_for_truncation = $pileup_str;
	
	if( @quals !=  $dep or $pileup_str =~ /\^\^/) {
		print STDERR join("\t",  ( $dep, $pileup_str, $qual_str  )), "\n\n";
		die;
	}
	
	if ( length ($pileup_str) == $dep) {
		@pileups = split "", $pileup_str;
		return @pileups;
	}
		
	#my $i = 0;
	#while ( $i < (length ($pileup_str) - 1 ) ){
	
	my $last_caret = "";
	
	while (  length($pileup_for_truncation) > 0 ){
		
		if ($debug) {
			print STDERR "len: ", length($pileup_for_truncation), "\n";#code
		}
				
		my $this_char = substr($pileup_for_truncation, 0, 1);
		my $next_char = substr($pileup_for_truncation, 1, 1);
		if ($this_char =~ /[ACGTNacgtn,.]/ and $next_char =~ /[^+-]/) { #not INDEL
			
			if ($debug) {
				print STDERR "1\n";
			}
			
			if ( $next_char eq "\$" ) { # end of a read
				if ($debug) {
					print STDERR  "2\n";
				}
				
				push @pileups, $last_caret . $this_char . $next_char;
				#$i+=2;
				$pileup_for_truncation = substr ( $pileup_for_truncation, 2 );
				$last_caret = "";
				next;
			}else{ # not end
				if ($debug) {
					print STDERR "3\n";
				}
				push @pileups, $last_caret . $this_char;
				#$i++;
				$pileup_for_truncation = substr ( $pileup_for_truncation, 1 );
				$last_caret = "";
				next;
			}
		}
		
		elsif( $this_char =~ /[ACGTNacgtn,.]/ and $next_char =~ /[+-]/){# INDEL
			#[+-][0-9]+[ACGTNacgtn]+
			if ($debug) {
				print STDERR "4\n";
			}
			if ( $pileup_for_truncation =~ /[ACGTNacgtn,.][+-](\d+)[ACGTNacgtn]/ ) {
				my $num = $1;
				my $digit_num = length($num);
				
				if ($debug) {
					print STDERR "num: $num, digit_num: $digit_num $this_char  $next_char \n\n";
				}
				
#				if ( $pileup_for_truncation =~ /([ACGTNacgtn,.][+-]\d{$digit_num,$digit_num}[ACGTNacgtn]{$num,$num})(.+)/ ) {
				if ( $pileup_for_truncation =~ /([ACGTNacgtn,.][+-]\d{$digit_num,$digit_num}[ACGTNacgtn]{$num,$num}\$?)(.+)/ ) {
					my ($first, $rest ) = ($1, $2);
					push @pileups, $last_caret . $first;
					$pileup_for_truncation = $rest;
					$last_caret = "";
					next;
				}
#				elsif( $pileup_for_truncation =~ /([ACGTNacgtn,.][+-]\d{$digit_num,$digit_num}[ACGTNacgtn]{$num,$num})$/ ) {
				elsif( $pileup_for_truncation =~ /([ACGTNacgtn,.][+-]\d{$digit_num,$digit_num}[ACGTNacgtn]{$num,$num}\$?)$/ ) {
					push @pileups, $last_caret . $pileup_for_truncation;
					last;				
				}
				
				else{
					print STDERR "4_else\n";
					print STDERR join("\t",  ("l178", $dep, $pileup_str, $qual_str, "trunc: " . $pileup_for_truncation  )), "\n\n";
					
					print STDERR join("\n", @pileups), "\n\n\n";
					
					die;
				}
			}
			else{
				print STDERR join("\t",  ("l185", $dep, $pileup_str, $qual_str, "trunc: " . $pileup_for_truncation  )), "\n\n";
				die;
			}
		}
		
		elsif( $this_char eq "*" ){
			if ($debug) {
				print STDERR "5\n";
			}
			push @pileups, $last_caret . $this_char;
			#$i++;
			$pileup_for_truncation = substr ( $pileup_for_truncation, 1 );
			$last_caret = "";
			next;
		}
		
		elsif( $this_char eq "^" ){
			if ($debug) {
				print STDERR "6\n";
			}
			$last_caret = substr ( $pileup_for_truncation, 0, 2 );
			$pileup_for_truncation = substr ( $pileup_for_truncation, 2 );
#			if ($pileup_for_truncation =~ /^\*/) {
#				print STDERR "carot cannot followed by asterisk\n\n";
#				print STDERR join("\t",  ( $dep, $pileup_str, $qual_str, "trunc: " . $pileup_for_truncation  )), "\n\n";
#				die;
#			}
			next;
		}
		elsif ($this_char =~ /[ACGTNacgtn,.]$/ ) {
			push @pileups, $last_caret . $pileup_for_truncation;
			last;
		}
		
		else{
			if ($debug) {
				print STDERR "7\n";
			}
			print STDERR join("\t",  ( "l219", $dep, $pileup_str, $qual_str, "trunc: " . $pileup_for_truncation  )), "\n\n";
			die;
		}
	}
	
	if ($print_results) {
		print STDERR  join("\n", @pileups), "\n\n";
		print STDERR  join("\n", @quals), "\n\n";
	}
	
	if ( @pileups !=  @quals){
		print STDERR "num not equal\n\n";
		print STDERR  join("\n", @pileups), "\n\n";
		print STDERR  join("\n", @quals), "\n\n";
		die;
	}
	return (@pileups);
}

#Kai_Module::determine_var($dep, \@seqs);

#^K	,	.	

sub determine_var{
	my ($dep, $ref) = @_;
	die unless (@{$ref} == $dep);
	
	my %records;
	
	for my $i(0..( $dep - 1 )){
		
		my $char = $ref->[$i];
		
		if ( $char =~ /^\^.[ACGTNacgtn.,]/ ) {
			$char = substr($char, 2);
		}
		elsif( $char =~ /\$$/ ){
			$char = substr($char, 0, length ($char)  - 1);
		}
		
		
		if ( $char eq  "*" ) {
			$records{ "DEL" } ++;
			next;
		}
		
		if ( length ($char)  == 1) {
			if ( $char =~ /[nN]/ ){
				$records{","}++;
			}else{
			    $records{ uc($char) } ++;
			}
			next;
		}
		
		
		
		
			
		if ( $char =~ /.\+\d+/ ) {
			if ( $char =~ /[,.]\+\d+/ ) {#INS
				$records{ "INS" } ++;
			}
			elsif( $char =~ /[ACGTNacgtn]\+\d+/ ){#SNP + INS
				$records{ "INS" } ++;	
			}
			else{
				print STDERR "l287\n";
				print STDERR join("\n", @{$ref}), "\n\n";
				die;
			}
		}		
			
		elsif($char =~ /.-\d+/){
			if ( $char =~ /[,.]-\d+/ ) {#DEL
				$records{ "DEL" } ++;
			}
			elsif( $char =~ /[ACGTNacgtn]-\d+/ ){#SNP + DEL
				$records{ "DEL" } ++;
			}
			else{
				print STDERR "l302\n";
				print STDERR join("\n", @{$ref}), "\n\n";
				die;	
			}
		}
		else{
			print STDERR "l308\n";
			print STDERR join("\n", @{$ref}), "\n\n";
			die;
		}
	}
##DEP=XX;TYPE=(SNP,INS,DEL);ALT=(ACTGID);PER=XX.XX,

	if ( check_hash_sum($dep, \%records) ){
		foreach my $k (sort keys %records){
			print STDERR join("\t", ($k, $records{$k})), "\n"
		}
		print STDERR "dep = $dep\n\n";
		die;
	}
	
	delete $records{"."};
	delete $records{","};
	
	my @keys = keys %records;
	if ( @keys == 0 ) {
		my $dep_item = "DEP=" . $dep;
		my $type_item = "TYPE=" . "NA";
		my $return_line = join(";", ($dep_item, $type_item ));
		return $return_line;
	}
	
	
	my @types;
	my @alts;
	my @pers;
	
	foreach my $k (sort keys %records){
		if ( $k =~/[ACTG]/ ) {
			push @types, "SNP";
			push @alts, $k;
			push @pers, sprintf("%.2f", 100 * $records{$k} / $dep );
		}
		elsif( $k eq "INS"){
			push @types, $k;
			push @alts, $k;
			push @pers, sprintf("%.2f", 100 * $records{$k} / $dep );
			
		}
		elsif( $k eq "DEL"){
			push @types, $k;
			push @alts, $k;
			push @pers, sprintf("%.2f", 100 * $records{$k} / $dep );
			
		}
		else{
			print STDERR $k, "\n";
			die;
		}
	}
	my $dep_item = "DEP=" . $dep;
	my $type_item = "TYPE=" . join(",", @types);
	my $alt_item = "ALT=" . join(",", @alts);
	my $per_item = "PER=" . join(",", @pers);
	my $return_line = join(";", ($dep_item, $type_item,$alt_item, $per_item));
	
	return $return_line;
	
}

#check_hash_sum($dep, \%records);
sub check_hash_sum{
	my ($dep, $ref) = @_;
	my $s = 0;
	foreach my $k(keys %{$ref}){
		$s+= $ref->{$k};
	}
	if ($s != $dep) {
		return 1;
	}
	else{
		return 0;
	}
}

#Kai_Module::read_list_recored_pos_and_list($input_list_file, \%pos_in_list_h, \@input_list, $chr_col );
sub read_list_recored_pos_and_list{
	my ($file, $pos_h, $list_a_ref, $chr_col) = @_;
	die unless (-e $file);
	open(IN , "$file") or die "cannot open $file";
	
	my $i = -1;
	
	while ( <IN> ) {
		$i++;
		chomp;
		$list_a_ref->[$i] = $_;
		next unless ($i>=1);
		my @a = split "\t";
		
		my ($chr, $s, $e ) = (@a[($chr_col - 1)..($chr_col + 1)]);
		$chr = simple_chr($chr);
		for my $pos( $s..$e){
			if (defined $pos_h->{$chr}->[$pos] ) {
				$pos_h->{$chr}->[$pos]  = $pos_h->{$chr}->[$pos]  . "," . $i;
			}else{
				$pos_h->{$chr}->[$pos]  = $i;
			}
			
		}
	}
	
	
	
	
	close IN;
}

#Kai_Module::cal_meth_level ( \@wmC_print, \%seqed_mC_h, \%seqed_dep_h, $i);
#Kai_Module::cal_meth_level ( \@mmC_print, \%seqed_per_h, \%covered_C_num_h, $i);
#Kai_Module::cal_meth_level ( \@fmC_print, \%seqed_isMeth_h, \%covered_C_num_h, $i);

sub cal_meth_level{
	my ($res_ref, $numerator_ref, $denominator_ref , $region_ind) = @_;
	my @types = ("CG", "CHG", "CHH", "C" );
	for my $i (0..$#types){
		my $t = $types[$i];
		if (defined $denominator_ref->{$region_ind}->{$t}) {
			#$ref_a->[$i] = $denominator_ref->{$region_ind}->{$t};
			my $x = 0;
			if (defined $numerator_ref->{$region_ind}->{$t} ){
				$x = $numerator_ref->{$region_ind}->{$t}
			}
			
			my $y = $denominator_ref->{$region_ind}->{$t};
			my $formula = "$x/$y";
			#Formula
			my $val = eval sprintf("%.2f", 100* $x/$y);
			$res_ref->[2*$i] = $formula;
			$res_ref->[2*$i + 1 ] = $val;
			die if $val > 100.01;
			
		}
		
	}
}

#	Kai_Module::cal_covered_per( \@covered_per_print, \@covered_C_num_print, \@C_num_print  );
sub cal_covered_per{
	my ( $res_ref, $numerator_ref, $denominator_ref ) = @_;
	my $last_index = scalar(@{$res_ref}) - 1;
	for my $i(0..$last_index){
		next if ( $denominator_ref->[$i] == 0);
		$res_ref->[$i] = eval sprintf ("%.2f", 100 * ($numerator_ref->[$i] ) / ( $denominator_ref->[$i] ) );
		die if $res_ref->[$i] > 100.01;
	}
}

#fill_C_num(\@C_num_print, \%C_num_h, $i);
sub fill_C_num{
	my ($ref_a, $ref_h, $region_ind) = @_;
	my @types = ("CG", "CHG", "CHH", "C" );

	for my $i (0..$#types){
		my $t = $types[$i];
		if (defined $ref_h->{$region_ind}->{$t}) {
			$ref_a->[$i] = $ref_h->{$region_ind}->{$t};
		}
		
	}
}


#Kai_Module::cal_win_num( \%num_win, \%chr_len, $win_size, $sliding_bp);

sub cal_win_num{
	my ($ref_num, $ref_len, $win_size_sub,  $sliding_bp_sub) = @_;
	
	foreach my $k (keys %{$ref_len}){
		my $len = $ref_len->{$k};
		my $n = int ( ( $len - 1 ) / $sliding_bp_sub ) + 1;
		$ref_num->{$k} = $n;
	}
}

#Kai_Module::hash_cumsum_SE ( \%cumsum_win,  \%num_win );

sub hash_cumsum_SE{
	my ($ref_res, $ref_in) = @_;
	my $sum = 0;
	foreach my $k(sort keys %{$ref_in}){
		my $s = $sum+1;
		my $e = $sum + $ref_in->{$k};
		$ref_res-> {$k} =  [ $s, $e];
		$sum += $ref_in->{$k};
	}
}

#my $ind_of_maximum = Kai_Module::get_ind_of_maximum(\@pers);
sub get_ind_of_maximum{
	my ($ref) = @_;
	my $last_index = scalar (@{$ref}) - 1;
	
	my $ind = 0;
	my $max = $ref->[$ind];
	for my $i(1..$last_index){
		if ( $ref->[$i] > $max ) {
			$max =  $ref->[$i];
			$ind = $i;
		}
	}
	return ($ind);
}


sub read_wig_hash{
	my ($wig_input, $hash_ref) = @_;
	print STDERR "read $wig_input...\t";
	open (IN, $wig_input) or die "cannot open $wig_input:$!";
	my $chr = 0;
	while (<IN>){
		chomp;
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = simple_chr( $1 );
			next;
		}
		
		my ($pos, $val) = split /\t/;
		$hash_ref->{$chr}->{$pos} = $val;
	}	
	close(IN);
	print STDERR "DONE\n";
}

sub read_wig_array{
	my ($wig_input, $hash_ref) = @_;
	print STDERR "read $wig_input...\t";
	open (IN, $wig_input) or die "cannot open $wig_input:$!";
	my $chr = 0;
	while (<IN>){
		chomp;
		next if (/^track/);
		if (/variableStep\s+chrom=(\w+)/){
			$chr = simple_chr( $1 );
			next;
		}
		
		my ($pos, $val) = split /\t/;
		$hash_ref->{$chr}->[$pos] = $val;
	}	
	close(IN);
	print STDERR "DONE\n";
}

#Feb 17, 2014
# for step1_extract_ChIP_read_num_pairs_using_featureCounts.pl
# input a list, output SAF file format for featureCounts
#read_list_and_get_sliding_interval ( $input, $type, $SAF_file, $bin_size, $sliding_size, $half_flanking_bin_num );

sub read_list_and_get_sliding_interval{
	my ( $input, $type, $output, $bin_size, $sliding_size, $half_flanking_bin_num ) = @_;
	open(IN, $input) or die "cannot open input";
	open(OUT, ">>$output") or die"cannot output";

	my $ID_index = 0;

	while (<IN>) {
		chomp;
		my @a = split "\t";
		next if ($a[0] =~ /(^coor|^#?chr$)/i );
		my ($chr, $start, $end);
		
	
	#	$ID_index++;
	
		if ($type eq "bed") {
			($chr, $start, $end) = @a[0..2];
			$chr = simple_chr($chr);
		}elsif($type eq "coor"){
			if ($a[0] =~ /(\S+):(\d+)-(\d+)/) {
				($chr, $start, $end) =  ($1, $2, $3);
				$chr = simple_chr($chr);
			}
		}else{
			die "wrong type\n\n";
		}
	
		my $mid_point = int ( ($start + $end)/2);
		my $first_start = $mid_point - int( ($bin_size - 1) /2 ) - $half_flanking_bin_num * $sliding_size  ;
	
		if($first_start <=0){
			print STDERR "$chr, $start, $end is close to start \n";
			next;
		}
		$ID_index++;

#	----------|----------
#	         + ++
		for my $bin_order ( (-1 * $half_flanking_bin_num)..$half_flanking_bin_num ){
			my $interval_start ;
				
			$interval_start = $mid_point - int( ($bin_size - 1) /2 ) +  $bin_order * $sliding_size ;
		
			my $interval_end   = $interval_start + $bin_size - 1 ; # exclude half read length
			if ($interval_start <= 0) {
				die $_;			
			}
		
			my $id =  $ID_index . ":" . $bin_order;
			print OUT join("\t", ($id, $chr, $interval_start, $interval_end, "+")), "\n";
		}	
	}

	close IN;
	close OUT;
}


#Feb 17, 2014
# for step1_extract_ChIP_read_num_pairs_using_featureCounts.pl
# input SAF file and bam dir

#featureCounts -T 8 -p -O -F SAF -a TAIR10_WinSize1000bp_sliding500bp_SAF.txt   -o  TAIR10_WinSize1000bp_sliding500bp_featureCounts_for_MBD7_ChIP_pO.txt  [bam]
#Kai_Module::run_featureCounts($SAF_file, $featureCounts_outfile, $feature_log_file, $CPU_thread, $bam_dir, $debug);

sub run_featureCounts{
	my ($SAF_file, $featureCounts_outfile, $feature_log_file, $CPU_thread, $bam_dir, $debug_featureCounts) = @_;
	unless ( -e $SAF_file){
		print STDERR "$SAF_file do NOT exists\n\n";
		die;
	}
	if ( -e $featureCounts_outfile) {
		print STDERR " $featureCounts_outfile exists \n\n";
		die;
	}
	unless ( -d $bam_dir) {
		print STDERR "$bam_dir do NOT exists\n\n";
	}
	
	opendir (DIR, $bam_dir) or die;
	my @bam_names = grep /\.bam$/, readdir DIR;
	close DIR;
	my @real_bams = map {File::Spec->catfile($bam_dir, $_)} @bam_names;
	
	my $bam_files = join(" ", @real_bams);
	my $cmd = "featureCounts -T $CPU_thread -p -O -F SAF -a $SAF_file -o $featureCounts_outfile $bam_files 2>> $feature_log_file";
	print STDERR $cmd, "\n\n";
	unless($debug_featureCounts){
		`$cmd`;
	}
}


#Gene_all_coordinate_SAF.txt
#"featureCounts" "-T" "3" "-p" "-O" "-F" "SAF" "-a" "Gene_all_coordinate_SAF.txt" "-o" "all_Gen
#Geneid  Chr     Start   End     Strand  Length  MBD7_5_My..bam
#0	1	 2	 3	 4	 5		6

#Kai_Module::extract_featureCounts_results($featureCounts_outfile, $output, $half_flanking_bin_num, $debug);

sub extract_featureCounts_results{
	my ( $featureCounts_outfile, $output, $half_flanking_bin_num,  $debug_extract_featureCounts ) = @_;
	die " $featureCounts_outfile Not Exists\n\n " if ( !$debug_extract_featureCounts and !(-e $featureCounts_outfile));
	
	if ( $debug_extract_featureCounts ) {
		print STDERR "debug : return from extract_featureCounts_results\n\n";
		return;
	}
	
	if (-e $output){
		print STDERR "output $output should not exists!!\n";
		die;
	}
	open(OUT, ">>$output") or die;
	open(IN, $featureCounts_outfile) or die;
	my $cmd_line = <IN>;
	my $head = <IN>;
	chomp $head;
	
	my @heads_a = split "\t", $head;
	my $names_start_index = 6;
	my @bam_names;# = @heads_a[$names_start_index..$#heads_a];
	
	#my ($volume, $directories, $infile_name) =       File::Spec->splitpath( $input );
	
	for my $i( $names_start_index..$#heads_a){
		my ($volume, $directories, $infile_name) =       File::Spec->splitpath( $heads_a[$i] );
		$infile_name =~ s/\.bam//;
		push @bam_names, $infile_name;
	}
	
	my %sums;
	my $num = 0;
	while (<IN>) {
		chomp;
		
		my @a = split "\t";
		
		#my $id =  $ID_index . ":" . $bin_order;
		my ( $ID_index , $bin_order ) = split ":", $a[0];
		$num = $ID_index;
		for my $i($names_start_index..$#a){
			$sums{($bam_names[($i - $names_start_index)])}->{$bin_order} += $a[$i];
		}
	}
	
#	print STDERR "file lines: $num\n";
	print OUT join("\t", ("Sample", (-1*$half_flanking_bin_num)..$half_flanking_bin_num)), "\n";
	
	foreach my $sample (sort keys %sums){
		print OUT $sample, "\t";
		
		for my $i ( (-1*$half_flanking_bin_num)..$half_flanking_bin_num) {
			my $tmp = 0;
			if (defined $sums{$sample}->{$i} ) {

				$tmp = eval sprintf ("%.2f",  ($sums{$sample}->{$i}) / $num );
				print OUT $tmp;
			}
			if ($i == $half_flanking_bin_num) {
				print OUT "\n";#code
			}else{
				print OUT "\t";
			}
			
			
		}
		
	}
	
	
	close IN;
	close OUT;
}

1;
__END__

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

if ($debug) {
	print STDERR , "\n";
}

http://www.perlmonks.org/?node_id=102347

First we get a namespace by declaring a package name. This helps ensure our module's
functions and variables remain separate from any script that uses it.

Use strict is a very good idea for modules to restrict the use of global variables.
See use strict warnings and diagnostics or die for more details.

We need to use the Exporter module to export our functions from the MyModule::
namespace into the main:: namespace to make them available to scripts that 'use' MyModule.

We pacify strict with the use vars declaration of some variables. We can use an 'our' declaration in 5.6+

We now set a $VERSION number and make Exporter part of MyModule using the @ISA. See perlboot for all the
gory details on what @ISA is or just use it as shown.

@EXPORT contains a list of functions that we export by default, in this case nothing. Generally the less
you export by default using @EXPORT the better. This avoids accidentally clashing with functions defined
in the script using the module. If a script wants a function let it ask.

@EXPORT_OK contains a list of functions that we export on demand so we export &func1 &func2 only if
specifically requested to. Use this in preference to just blindly exporting functions via @EXPORT.
You can also export variables like $CONFIG provided they are globals not lexicals scoped with my
(read declare them with our or use vars).

%EXPORT_TAGS. For convenience we define two sets of export tags. The ':DEFAULT' tag exports only &func1;
the ':Both' tag exports both &func1 &func2. This hash stores labels pointing to array references.
In this case the arrays are anonymous.

We need the 1; at the end because when a module loads Perl checks to see that the module returns
a true value to ensure it loaded OK.
You could put any true value at the end (see Code::Police) but 1 is the convention.


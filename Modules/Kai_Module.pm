package Kai_Module;


use strict;
use warnings;
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
			$records{ uc($char) } ++;
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

1;
__END__

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


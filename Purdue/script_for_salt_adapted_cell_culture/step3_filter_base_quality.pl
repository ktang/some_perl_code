#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;


my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

my $min_base_qual = 20;


my $usage = "$0 \n <input> <output> <sample_num> <Phred_33_64>\n\n";
die $usage unless(@ARGV == 4);

my $input = shift or die;
my $output = shift or die;
my $sample_num  = shift or die;
my $phred_33_64 = shift or die;
 
die "Phred_33_64" unless ($phred_33_64 == 33 or $phred_33_64 == 64);

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $head = <IN>;
print OUT $head;

while (<IN>) {
	chomp;
	my @a = split "\t";
	print OUT join("\t", (@a[0..2])), "\t";
	for my $n ( 1..$sample_num ){
		
		my ($old_dep, $old_seq, $old_qual) = @a[ ( 3 * $n )..( 3*$n + 2)];
		
		if ( $old_dep == 0) {
			print OUT join("\t", ($old_dep, $old_seq, $old_qual  )) ;
		}
		else{
			my @seqs = Kai_Module::split_pileup( $old_dep, $old_seq, $old_qual );
			my @quals = split "", $old_qual;
			if (@seqs != @quals) {
				print STDERR "seqs: ", @seqs, "\n\n";
				print STDERR "quals:", @quals, "\n\n";
				die;
			}
		
			my $new_seqs =  "";
			my $new_quals = "";
			my $new_dep   = 0;
		
			for my $i (0..$#quals){
				my $this_seq  = $seqs[$i];
				my $this_qual = $quals[$i];
				die unless ( length($this_qual) == 1);
				if ( (ord($this_qual) - $phred_33_64 ) >=  $min_base_qual) {
					$new_dep ++;
					$new_seqs  .= $this_seq;
					$new_quals .= $this_qual;
				}
			}
			
			if ($new_dep == 0) {
				print OUT join("\t",  (0, "*", "*"));
			}
			else{
				my @tmp_q = split "", $new_quals;
				my @tmp_p = Kai_Module::split_pileup ($new_dep, $new_seqs, $new_quals  );
				
				if ( @tmp_q != $new_dep or @tmp_p != $new_dep) {
					print STDERR "l89 seqs: ", @seqs, "\n\n";
					print STDERR "quals:", @quals, "\n\n";
					die;
				}
				
				
				print OUT join("\t", ($new_dep, $new_seqs, $new_quals ) ) ;
			}
		}
		if ($n == $sample_num) {
			print OUT  "\n";
		}else{
			print OUT "\t";	
		}
	}
	
}



close(IN);
close(OUT);

exit;

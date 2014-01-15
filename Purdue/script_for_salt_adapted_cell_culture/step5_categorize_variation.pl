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

#my $min_base_qual = 20;

my $min_dep = 8;


my $usage = "$0 \n <input> <output> <sample_num> \n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $output = shift or die;
my $sample_num  = shift or die;
#my $phred_33_64 = shift or die;
 
#die "Phred_33_64" unless ($phred_33_64 == 33 or $phred_33_64 == 64);

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my $head = <IN>;
#print OUT $head;

chomp $head;
my $head_out = deal_head($head, $sample_num);
print  OUT $head_out, "\n";

#DEP=XX;TYPE=(SNP,INS,DEL);ALT=(ACTGID);PER=XX.XX
while (<IN>) {
	chomp;
	my @a = split "\t";
#	print OUT join("\t", (@a[0..2])), "\t";
#Kai_Module
	
	print OUT join ("\t", @a[0..2]), "\t";
	
	for my $n ( 1..$sample_num ){
		my ( $dep, $seq, $qual ) = @a[ ( 3 * $n )..( 3*$n + 2)];
		
		if ( $dep == 0 ) {
			print OUT "DEP=0";
		}
		else{
			my @seqs = Kai_Module::split_pileup( $dep, $seq, $qual );
			die unless ( @seqs == $dep);
			my $var = Kai_Module::determine_var($dep, \@seqs);
			print OUT $var;
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

#my $head_out = deal_head($head);
sub deal_head{
	my ($line, $num) = @_;
	my @a = split "\t", $line;
	my @tmp = ();
	
	for (my $i = 3; $i<= 3* $num; $i+= 3){
		my $item = $a[$i];
		if ( $item =~ /[a-zA-Z]+_(\S+)/) {
			my $label = $1;
			push @tmp, $label;
		}
		
	}
	
	my $out_line = join("\t", (@a[0..2], @tmp));
	return $out_line;
	
}
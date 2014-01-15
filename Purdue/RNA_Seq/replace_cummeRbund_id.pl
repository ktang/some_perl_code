#!/usr/bin/perl -w

# purpose :
# the output of cummeRbund looks like this:

#gene_id	gene_id.1	sample_1	sample_2	status	value_1	value_2	log2_fold_change	test_stat	p_value	q_value	significant
#XLOC_000035	XLOC_000035	FQ_WT	X58_2_ibm2	OK	7.88784	0.987371	-2.99797	3.70059	0.000215096	0.0099636	yes
#XLOC_000068	XLOC_000068	FQ_WT	X58_2_ibm2	OK	12.6479	112.855	3.1575	-4.18719	2.82426e-05	0.00207208	yes

# replace XLOC_000035 with gene AGI

#gene_id	has_more_isoform

use strict;
use File::Spec;

my $db = "/Volumes/My_Book/20121228_RNA_seq/clout/diff_out/isoforms.fpkm_tracking";
die unless (-e $db);

my $debug = 0;

my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

my ( %ids, %more_isoforms);

read_db($db, \%ids, \%more_isoforms);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">>$output") or die "cannot open $output: $!";

my $head = <IN>;
chomp $head;
my @a_h = split "\t", $head;

$a_h[1] = "has_more_isoforms";
print OUT join("\t", @a_h[0..$#a_h-1]), "\n";

while(<IN>){
	chomp;
	my @a = split "\t";
	die $_ unless (defined $ids {$a[0]});
	
	$a[1] = "no";

	if(defined $more_isoforms {$a[0]}){
		$a[1] = "yes";
	}
	
	$a[0] = $ids {$a[0]} ;
	
	print OUT join("\t", @a[0..$#a-1]), "\n";
	
}

close(IN);
close(OUT);

exit;

#read_db($db, \%ids, \%more_isoforms);
sub read_db{
	my ($file, $id_ref, $isoform_ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die;
	my $head = <IN>;
	while (<IN>){
		chomp;
		my @a = split "\t";
		my $xloc_id = $a[3];
		my $gene_id = $a[2];
		if (defined $id_ref -> { $xloc_id }){
			my $t = $id_ref -> { $xloc_id };
			if( $gene_id ne $t){
				$isoform_ref-> {$xloc_id } = 1;
			}
		}else{
			$id_ref -> { $xloc_id } = $gene_id;
		}
	}
	close(IN);
	
}

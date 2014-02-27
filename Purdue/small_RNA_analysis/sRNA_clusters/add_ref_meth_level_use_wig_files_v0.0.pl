#!/usr/bin/perl -w
use strict;
use File::Spec;

# add_ref_meth_level_use_wig_files_v0.0.pl
# purpose: use the wig file ( like /Volumes/Macintosh_HD_2/Jacobsen_Cell_2013/wig/all_C/GSM980986_WT_rep2_mC.wig )
# to calculate methylation level like mmC

my $debug = 0;
if($debug){
	print STDERR "debug = 1\n\n";
}

my $usage = "$0 \n <input_bed_like_file> <wig> <label> <output>\n\n";
die $usage unless(@ARGV == 4);


my $bed_file = shift or die;
my $wig_file = shift or die;
my $label = shift or die;
my $output = shift or die;

die unless (-e $bed_file) ;
die unless (-e $wig_file);
die if(-e $output);

if ($debug) {
	print join("\n", ($bed_file, $wig_file, $output)), "\n\n";
	exit;
}



my @bed_list;
my %pos;


my $last_one_index = read_bed_file($bed_file, \@bed_list);
record_region_pos(\@bed_list, \%pos);

my @bed_wig_val;

push_wig_val($wig_file, \%pos, \@bed_wig_val);

output_file($output, \@bed_list, \@bed_wig_val, $label);


exit;

# output_file($output, \@bed_list, \@bed_wig_val);
sub output_file{
	my ($file, $list_ref, $val_ref, $label_ref) = @_;
	die if(-e $file);
	
	my $last_index = scalar(@{$list_ref}) - 1;

	
	open(OUT, ">$output") or die;
	my $raw_head = $list_ref->[0];
	print OUT join("\t", ( $raw_head, "wig_formula_" . $label_ref, "wig_mmCXX_" . $label_ref )), "\n";
	foreach my $i(1..$last_index){
		my $raw_line = $list_ref->[$i];
		my ($formula, $level) = ("NA", "NA");
		if ( defined $val_ref->[$i]) {
			my @tmp = @{$val_ref->[$i]};
			($formula, $level) = cal_level(\@tmp);
		}
		print OUT join("\t", ($raw_line, $formula, $level)), "\n";
	}
	
	close OUT;
}

sub cal_level{
	my ($ref) = @_;
	my $num = scalar(@{$ref});
	my $sum = 0;
	foreach my $i (@{$ref}){
		$sum+= $i;
	}
	$sum*=100;
	return ( "$sum/$num", sprintf("%.3f", $sum/$num));
}

sub push_wig_val{
	my ($wig_sub, $pos_ref, $bed_wig_val_ref) = @_;
	die unless (-e $wig_sub);
	open(IN, $wig_sub) or die "wig";
	my $chr;
	while (<IN>) {
		next if(/^track/);
		if (/chrom=(\S+)/) {
			$chr = lc($1);
			next;
		}
		chomp;
		my ($pos, $val) = split "\t";
		$val = abs($val);
		if (defined $pos_ref->{$chr}->[$pos]) {
			my $index = $pos_ref->{$chr}->[$pos];
			push @{$bed_wig_val_ref->[$index]}, $val;
		}
 	}
	close IN;
}

sub read_bed_file{
	my ($file, $ref) = @_;
	die unless (-e $file);
	
	open(IN, $file) or die;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		$ref->[$i] = $_;
	}
	close(IN);
	return $i;
}

sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ($chr, $start, $end ) = @a[0..2];
		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i; #record index
		}
	}
}
#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_GFF> <output_exon_GFF>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

#open(IN, $input) or die "cannot open $input: $!";
#die if(-e $output);
#open(OUT, ">$output") or die "cannot open $output: $!";
#close(IN);
#close(OUT);


my (%strands, %gff_records,  %intron_records);

read_GFF ($input, \%strands, \%gff_records);

extract_intron(\%strands, \%gff_records, \%intron_records);

output_intron($output, \%intron_records, \%strands);


exit;

#Chr1    TAIR10  exon    3631    3913    .       +       .       Parent=AT1G01010.1


#read_GFF ($input, \%strands, \%gff_records);
sub read_GFF{
	my ($file, $strand_ref, $gff_records_ref) = @_;
	die unless (-e $file);
	open(IN, $file) or die "cannot open $file: $!";
	
	while(<IN>){
		next unless /exon/;
		chomp;
		my @a = split "\t";
		next unless ($a[2] eq "exon");
		die $_ unless ($a[2] eq "exon");
		my $temp;
		my $trascript;
		($temp, $trascript)  = split "=", $a[-1];
		die $_ unless ($trascript =~ /AT.G\d\d\d\d\d\.\d/);
		
		$strand_ref->{ $trascript } = join("_", ( lc($a[0]), $a[6]));
		my $interval = join("_", @a[3..4]);		
		push @{ $gff_records_ref ->{ $trascript }}, $interval;		
	}
	close IN;
}

#get_exon($exon1,  $exon2 );
sub get_exon{
	my ($exon1,  $exon2 ) = @_ ;
	
	my ($s1, $e1) = split "_", $exon1 ;
	my ($s2, $e2) = split "_", $exon2 ;
	my $intron_inverval = join("_", ($e1+1, $s2 - 1));
	return $intron_inverval;
}

#extract_intron (\%strands, \%gff_records, \%intron_records);
sub extract_intron{
	my ( $strand_ref, $gff_records_ref,  $intron_records_ref) = @_;
	
	foreach my $trascript ( keys %{$gff_records_ref} ){
		my @exon_array = @{ $gff_records_ref ->{ $trascript } };
		
		next if ( @exon_array == 1 );
		die $trascript unless (defined $strand_ref->{ $trascript });
		
		my $label = $strand_ref->{ $trascript };
		my ($chr, $strand ) = split "_", $label;
		
		#-->	-->	-->
		if( $strand eq "+" ){
			for my $i(0..($#exon_array - 1)){
				my $exon_left = $exon_array[$i];
				my $exon_right = $exon_array[$i + 1];
				my $intron_inverval = get_exon($exon_left,  $exon_right );
				push @{ $intron_records_ref ->{ $trascript }}, $intron_inverval;
			}
		}
		#2	1	0
		#<--	<--	<--
		elsif( $strand eq "-" ){
			for my $i(0..($#exon_array - 1)){
				my $exon_right = $exon_array[$i];
				my $exon_left = $exon_array[$i + 1];
				my $intron_inverval = get_exon($exon_left,  $exon_right );
				push @{ $intron_records_ref ->{ $trascript }}, $intron_inverval;
			}
		}else{
			print STDERR $label;
			die $trascript;
		}
		
	}
	
}

#output_intron($output, \%intron_records);
sub output_intron{
	my ($file, $intron_records_ref, $strand_ref ) = @_;
	
	die if (-e $file);
	open(OUT, ">>$file") or die;
	
	foreach my $trascript ( sort keys %{$intron_records_ref} ){
		die $trascript unless (defined $strand_ref->{ $trascript });
		
		my $label = $strand_ref->{ $trascript };
		my ($chr, $strand ) = split "_", $label;
		
		my @intron_array = @{ $intron_records_ref ->{ $trascript } };
		
		for my $i (0..$#intron_array){
			my $interval = $intron_array[$i];
			my ($start, $end) = split "_", $interval;
			print OUT join("\t", ($chr, "TAIR10", "intron", $start, $end, ".", $strand, ".", "Parent=". $trascript )), "\n";
		}
		
	}
	
	close(OUT);
	
}
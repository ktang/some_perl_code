#!/usr/bin/perl -w
# Copyright (C)  2010-
# Program:			extract_marker_pos_in_consensus.pl
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		Jan 18,2011 @Purdue
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	Jan 26,2011
# Description:	
#**************************
# Version: 1.0	read SNP file and store in a hash, then read the consensus file to get
#				the information in the SNP position.
# Version: 1.1	add two col in the output, one is main_SNP/(ref + main_SNP) 
#				the last is main_SNP/total count
#**************************
# e-mail:tangkai.ustc@gmail.com

=head2
	file format:
	SNP marker file
	<strain_name>.SNPs.TAIRX.txt
              <Sample>
              <Chromosome>
              <Position>
              <Reference base>
              <Substitution base>
              <Quality>
              <# of nonrepetitive reads supporting substituion>
              <concordance>
              <Avg. # alignments of overlapping reads>
              
     consensus file
     Chr	pos	ref_base	depth	A	C	G	T
     
     output file
     Chr	pos	ref_base	ref_num	Ler_base	Ler_num	A	C	G	T	main_SNP/(ref + main_SNP)	main_SNP/total_count
=cut

use strict;
use Getopt::Long;



my $version="1.1";

my $marker_file = ""; # SNP marker file
my $consensus_file = "";
my $output_file = "";

my %CMD;
GetCom();

print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #ÔËÐÐ¿ªÊ¼Ê±¼ä
print STDERR "Now = $Time_Start\n\n";


### Parse snp marker file---------------------------------------------------------------------
my %pos_hash = ();
my %ref_h = ();
my %Ler_h = ();
open SNP,$marker_file or die "cannot open $marker_file:$!";

while (<SNP>){
	my @pts = split "\t";
	$pos_hash{"Chr".$pts[1]."#".$pts[2]} = 1;
	$ref_h{"Chr".$pts[1]."#".$pts[2]} = $pts[3];
	$Ler_h{"Chr".$pts[1]."#".$pts[2]} = $pts[4];
}


### Parse consensus file-----------------------------------------------------------------------

open CON, $consensus_file or die "Cannot open $consensus_file:$!";
open (OUT, ">$output_file")
	or die "Cannot open $output_file:$!";

while(<CON>){
	chomp;
	my @pts = split "\t";
	if (defined($pos_hash{$pts[0]."#".$pts[1]})){
		my %base = ();
		$base{"A"} = $pts[4];
		$base{"C"} = $pts[5];
		$base{"G"} = $pts[6];
		$base{"T"} = $pts[7];
		my $ref = $ref_h{$pts[0]."#".$pts[1]};
		my $Ler = $Ler_h{$pts[0]."#".$pts[1]};
		my $ref_num = $base{$ref};
		my $Ler_num = $base{$Ler};

##########change in v1.1
		my $per_1 = -0.01;
		my $per_2 = -0.01;
		if ($ref_num + $Ler_num != 0){
			$per_1 = $Ler_num/($ref_num + $Ler_num);
			$per_2 = $Ler_num/$pts[3];
		}
		
		else{
			print STDERR $_, "\n";
			next;
		}
##########change in v1.1		
		
		my $out_line = join("\t",$pts[0],$pts[1],$ref,$ref_num,$Ler,$Ler_num,$pts[4],$pts[5],$pts[6],$pts[7],$per_1,$per_2);
		if ($ref_num + $Ler_num != $pts[3]){
			print $out_line,"\n";
		}
		
		print OUT $out_line,"\n";
	}
}

sub_end_program();
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime 
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}


############################################################################################################
######################                  GetCom()
############################################################################################################
sub GetCom{
	my @usage = ("\nUsage: $0
	
	Mandatory:
	--marker	STRING	SNP marker file format:<Sample>	<Chromosome>	<Position>	<Reference base>	<Substitution base>
	--consensus	STRING	consensus file format:Chr	pos	ref_base	depth	A	C	G	T
	--output	STRING	output file name
	note: out put inconsistent postion in STDOUT	
	writtern by Kai Tang.
	Purdue University, West Lafayette, 2011.
	
	\n");
	
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD,"marker=s","consensus=s","output=s");
	
	die("Please specify snp marker file\n") unless defined($CMD{marker});
	die("Please specify consensus file\n") unless defined($CMD{consensus});
	die("Please specify output file\n") unless defined($CMD{output});
	
	$marker_file = $CMD{marker}; 
	$consensus_file = $CMD{consensus};
	$output_file = $CMD{output};
}

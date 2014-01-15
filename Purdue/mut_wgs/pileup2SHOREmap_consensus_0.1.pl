#!/usr/bin/perl -w
# Copyright (C)  2010-
# Program:			pileup2SHOREmap_consensus
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		Dec 15,2010
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	Jan 11, 2011
# Description:		convert pileup to SHOREmap input consensus
#**************************
# Version: 0.0
# Version: 0.1 focus on Chr1-5,ignore ChrC and ChrM.
#			yes|no indicate whether output header
#**************************
# e-mail:tangkai.ustc@gmail.com

#output format
#Chr	pos	ref	depth	A	C	G	T

my $version="0.1";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time()));
print STDERR "Now = $Time_Start\n\n";




my $usage = "$0 <pileup_input> <yes|no> stdout";

die $usage unless (@ARGV == 2);

my $pileup_file = $ARGV[0];
my $header = $ARGV[1];
die $usage unless ($header =~ /yes|no/);

if ($header eq "yes"){
	print join("\t",("Chr","position","ref_base","depth","A","C","G","T")),"\n";
}

open (IN, $pileup_file)
	or die "cannot open $pileup_file:$!";
	
while(<IN>){
	my @pts = split /\t/;
	
	
	if ($pts[0] eq "ChrC" or $pts[0] eq "ChrM"){
		next;
	}
	
	#################
	
	#################
	
	my %nucleotides; # a hash store the number of ATCG
#	my $num_ref = ( $pts[4] =~ tr/.,/.,/);
#	$nucleotides{$pts[2]} = $num_ref;
	$nucleotides{$pts[2]} = ($pts[4] =~ tr/.,/.,/);
	
	unless (defined $nucleotides{"A"}){
		$nucleotides{"A"} = ($pts[4] =~ tr/Aa/Aa/);
	}
	
	unless (defined $nucleotides{"C"}){
		$nucleotides{"C"} = ($pts[4] =~ tr/Cc/Cc/);
	}
	
	unless (defined $nucleotides{"G"}){
		$nucleotides{"G"} = ($pts[4] =~ tr/Gg/Gg/);
	}
	
	unless (defined $nucleotides{"T"}){
		$nucleotides{"T"} = ($pts[4] =~ tr/Tt/Tt/);
	}	
	
#	my @acgt = ("A","C","G","T");

#	foreach my $nt (@acgt){
#		if(defined $nucleotides{$nt}){
#			next;
#		}
#		else{
#			
#			$nucleotides{$nt} = ($pts[4] =~ tr/$nt/$nt/i)
#		}
#	}

	my $total = $nucleotides{"A"} + $nucleotides{"C"} + $nucleotides{"G"} + $nucleotides{"T"};
	
	if ($total != $pts[3]){
		print STDERR "$_",join ("\t",($nucleotides{"A"},$nucleotides{"C"},$nucleotides{"G"},$nucleotides{"T"})) , "\n";
		sub_end_program();
		die;
	}
	
	my @out_pts = ($pts[0],$pts[1],$pts[2],$pts[3],$nucleotides{"A"},$nucleotides{"C"},$nucleotides{"G"},$nucleotides{"T"});
	print join("\t",@out_pts),"\n";
}

sub_end_program();

############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);
}
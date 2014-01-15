#!/usr/bin/perl -w
# create wig file to display methlated positions
use strict;
use Math::CDF qw(:all);

my $label = shift @ARGV;
my $conver_error = 0.0275;
my $sig_level = 0.05;
#my $pos_out = "flowcell41_lane8_pairs_meth_plus.wig";
#my $neg_out = "flowcell41_lane8_pairs_meth_minus.wig";
#open(PO, ">$pos_out") or die "Can't open $pos_out: $!";
#open(NO, ">$neg_out") or die "Can't open $neg_out: $!";
print 'track type=wiggle_0 name="', $label, '" description="', $label, ' methylation" viewLimits=-1.0:1.0 color=255,102,0', "\n";

#print NO 'track type=wiggle_0 name="Col0_minus" description="Col0 methlation on minus strand" viewLimits=0.0:2.0 color=255,102,0', "\n";

#my %vals = ("CG"=>2, "CHG"=>1.5, "CHH"=>1);
my %meth; #chr->pos=strand
while(<>){
	chomp;
	my ($chr, $pos, $pos2, $total, $methy_level, $strand) = split /\t/;
	next if($chr eq "chrM" || $chr eq "chrC" || $chr eq "pUC19");
	next if($total < 2);
	my $numC = $methy_level * $total;
	my $p = pbinom($numC, $total, (1-$conver_error));
	if($p > $sig_level){ # methylated
		if($strand eq "+"){
		    $meth{$chr}->{$pos+1} = $methy_level;
		}else{
			$meth{$chr}->{$pos+1} = -1*$methy_level;
		}
	}
}

foreach my $chr(sort keys %meth){
	print "variableStep chrom=$chr span=1\n";
	my %data = %{$meth{$chr}};
	foreach my $p(sort {$a<=>$b} keys %data){
		print $p, "\t", $data{$p}, "\n";
		
		#if($data{$p} > 0){
		#	print $p, "\t", 1, "\n";
	    #}else{
		#	print $p, "\t", -1, "\n";
		#}
	}
}


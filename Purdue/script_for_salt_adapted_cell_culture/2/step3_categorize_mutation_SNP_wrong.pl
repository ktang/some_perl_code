#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );

# input and output is std
# cat *vcf | perl $0 > output.txt

use strict;
#less 6sample_mpileup_in_5000loci.txt | cut -f1-4,5,7,8,10,11,13,14,16,17,19,20 | less | grep -P '[\+-]' | less
# less 6sample_mpileup_in_5000loci.txt | cut -f1-4,5,7,8,10,11,13,14,16,17,19,20 | less | grep -v -P '[\+-]' | less

#chr1    3       C       3       ..^I.+7TAAACCC  1       .       1       .       0       *       1       .       0       *
#chr1    4       T       3       ...     1       .       2       .^IG+7AAACCCT   0       *       1       A+9CTAAACCCT    0       *
#0	 1	 2	 3	 4	 5	 6	 7		8	 9	 10	 11	 12		 13	 14
#                            WT0             WT150             ddc0                 ddc150          nrpe0                    nrpe150
my %records;

die "incompeleteed\n\n";
my %refs;

while (<>){
	chomp;
	my @a = split "\t";
	my $indel = 0;
	my ( $chr, $pos, $ref ) = @a[0..2];
	
	my @b;
	
	for  (my $i = 3; $i<=13; $i+=2){
		push @b, find_SNP($chr, $pos, $ref , $a[$i] , $a[$i + 1] );
	}
	
	print join("\t", ( $chr, $pos, $ref , @b )), "\n";
}

exit;

sub find_SNP{
	my ($chr_sub, $pos_sub, $ref_sub, $dep, $str) = @_;
	
	if( $dep == 0 ){
		return "NA";
	}
	
	my $dot = ( $str =~ tr/\./\./ );
	my $comma = (  $str =~ tr/,/,/   );
	my $A_num = ( $str =~ tr/Aa/Aa/  );
	my $C_num = ( $str =~ tr/Cc/Cc/  );
	my $G_num = ( $str =~ tr/Gg/Gg/  );
	my $T_num = ( $str =~ tr/Tt/Tt/  );
	
	my ($m1, $m2, $m3, $m4) = max_sort($A_num, $C_num, $G_num, $T_num);
	
	if(  $dot + $comma + $A_num + $C_num + $T_num + $G_num != $dep or $m2 > 0){
		print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $str, $dot , $comma , $A_num , $C_num , $T_num , $G_num )), "\n", 
	}
	my $max = $m1;
	my $hit_num = 0;
	my $snp_base = "N";
	if( $max == $A_num) {$hit_num ++; $snp_base = "A"; }
	if( $max == $C_num) {$hit_num ++; $snp_base = "C"; }
	if( $max == $G_num) {$hit_num ++; $snp_base = "G"; }
	if( $max == $T_num) {$hit_num ++; $snp_base = "T"; }
	
	if  ($m1 ==  $m2 and $m1 > 0) {
		print STDERR join("\t", ($max, $A_num, $C_num, $G_num, $T_num, $chr_sub, $pos_sub, $ref_sub, $dep, $str)), "\n" ;#if($max >= 3);
		$snp_base = "N";
	}
	
	my $per = sprintf("%.2f", $max /$dep *100);
	return ("$ref_sub=>$snp_base,$per");
}

sub max_sort{
	my @a = @_;
	my @b = sort {$b <=> $a} @a;
	return @b;
}

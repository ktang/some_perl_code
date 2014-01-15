#!/usr/bin/perl -w

#step3_categorize_mutation_SNP_v0.2.pl
# modified from 
#/Users/tang58/misc/Zhu_Xiaohong/downstream/src/step3_categorize_mutation_SNP_v0.1.pl

# Apr 17, 2013
#purpose I found combined 513 loci showing variation.
# the program is to categorize them.

#      1 chr
#      2 pos
#      3 ref
#      4 dep_WT_0
#      5 seq_WT_0
#      6 qual_WT_0
#      7 dep_WT_150
#      8 seq_WT_150
#      9 qual_WT_150


# v0.1
# foreach locus categorize it

#v0.0 uncompleted Mar10, 2013

#input
#~/misc/Zhu_Xiaohong/downstream/dep5_13:19:47_N=586$
#less 6sample_mpileup_in_5000loci_filter_allDep5_3983loci.txt 

#chr1    88937   G       14      .,,.,,..,,.,.,  DHJHJJIIJJDJCB  15      A,a,,.aAAa,aAA. 7H1JHH2824>6:3C 15      AAAAAaAAaAaaAa^Ka       887=61.5134755! 5       AAAaa   7325!   15      AaAaAAAAaAaaaAA 72857633231376; 11      aAaa      3 AaAAAaA     5731618:2=6
# 0	  1	 2	 3		4		5; 6-8; 	9-11; 		12-14; 		15-17; 		19-20
#					WT0		    WT_150	 ddc_0		ddc_150		nrpe_0		nrpe150
					
use strict;

my $debug = 1;
my %records;

my %refs;

my $usage = "$0 \n <input_db_file> STDOUT\n\n";
die $usage unless (@ARGV == 1);
my $input = shift or die "input";
die unless (-e $input);

open (IN, $input) or die;

my $h = <IN>;
chomp $h;
my @tmp = split "\t", $h;

print join( "\t",  ( @tmp[0..2],
		    "mut_WT_0",		"per_WT_0",
		    "mut_WT_150",	"per_WT_150",
		    "mut_ddc_0",	"per_ddc_0",
		    "mut_ddc_150",	"per_ddc_150",
		    "mut_nrpe1_0",	"per_nrpe1_0",
		    "mut_nrpe1_150",	"per_nrpe1_150",
		     @tmp[3..$#tmp]
		    )), "\n";

while (<IN>){
	chomp;
	my @a = split "\t";
#	my $indel = 0;
	my ( $chr, $pos, $ref ) = @a[0..2];
	
	my @b = ();
	
	for  (my $i = 3; $i<=18; $i+=3){
	#	push @b, $a[$i];
		push @b, categorize($chr, $pos, $ref , $a[$i] , $a[$i + 1] , $i);
	}
	
	# print join("\t",  ( $chr, $pos, $ref , @b , @a[4,7,10,13,16,19]) ), "\n";
	
	print join("\t",  (@a[0..2], @b, @a[3..$#a ] )), "\n"
	
#	print join("\t", ( $chr, $pos, $ref , @b )), "\n";
}

close IN;

exit;

sub categorize{
	my ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub) = @_;
	my $str_modified = $str;
	
	if ($dep == 0){
		return ("NA", "NA");
	}
	
#	my @chars = ();
	
	my $arrow_num  = ( $str_modified =~ tr/\^/\^/ );
	if($arrow_num > 0){
	#	remove_char_and_one_after(\$str_modified, "\^", $arrow_num);
		#if($debug){ 		#	print STDERR "arrow_l51: $arrow_num, $chr_sub, $pos_sub\n";}
		remove_arrow( \$str_modified,  $arrow_num );
		
	}
	
#	my $dollar_num = ( $str_modified =~ tr/\$/\$/ );
#	if ($dollar_num > 0){
		#remove_char_and_one_after(\$str_modified, "\$", $dollar_num);
#		if($debug){#	print STDERR "dollar_num_l62: $dollar_num, $chr_sub, $pos_sub \n";}
#		remove_dollar ( \$str_modified,  $dollar_num  );
#	}
	
	my $plus  = ( $str_modified =~ tr/\+/\+/ );
	if($plus > 0){
		remove_insersion_chars(\$str_modified,  $plus, $chr_sub, $pos_sub, $i_sub);
	}
	
	my $minus = ( $str_modified =~ tr/-/-/ );
	if ( $minus > 0 ){
		remove_deletion_chars( \$str_modified, $minus, $chr_sub, $pos_sub, $i_sub);
	}
	
	
	my $dot_num   = ( $str_modified =~ tr/\./\./  );
	my $comma_num = (  $str_modified =~ tr/,/,/   );
	my $A_num = ( $str_modified =~ tr/Aa/Aa/  );
	my $C_num = ( $str_modified =~ tr/Cc/Cc/  );
	my $G_num = ( $str_modified =~ tr/Gg/Gg/  );
	my $T_num = ( $str_modified =~ tr/Tt/Tt/  );
	
	my $N_num = ( $str =~ tr/Nn/Nn/  );
	my $star_num = ( $str =~ tr/\*/\*/  );
	
	my $sum =  $dot_num + $comma_num + $A_num + $C_num + $T_num + $G_num  + $plus + $minus +$N_num + $star_num ;
	
#	if( $dep !=  $sum and $sum - $dep != 1){ # $arrow_num + $dollar_num
		#print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $str, $str_modified  ,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus,  $arrow_num , $dollar_num  )), "\n", 
	#	 print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $sum, $str, $str_modified  ,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus , $N_num , $star_num )), "\n", 
#	}
	
	my ($m1, $m2, $m3, $m4) = max_sort($A_num, $C_num, $G_num, $T_num, $plus, $minus + $star_num );
	
#	if(  $dot + $comma + $A_num + $C_num + $T_num + $G_num + $N_num + $star_num != $dep or $m2 > 0){
#		print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $str, $dot , $comma , $A_num , $C_num , $T_num , $G_num , $N_num , $star_num)), "\n", 
#	}

	if ($m2 == $m1 and $m1 > 1 ){
	#	print join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $sum, $str, $str_modified  ,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus , $N_num , $star_num )), "\n" 
	}
	
	my $max = $m1;
	my $hit_num = 0;
	my $snp_base = "N";
	
	my $per = sprintf( "%.2f" , $max / $dep * 100);
	
	if   ( $max == $A_num 		 ) {$hit_num ++; $snp_base = "A"; }
	elsif( $max == $C_num 		 ) {$hit_num ++; $snp_base = "C"; }
	elsif( $max == $G_num 		 ) {$hit_num ++; $snp_base = "G"; }
	elsif( $max == $T_num 		 ) {$hit_num ++; $snp_base = "T"; }
	elsif( $max == $plus  		 ) {$hit_num ++; $snp_base = "INS"; }
	elsif( $max == $minus + $star_num) {$hit_num ++; $snp_base = "DEL"; }
	else{
		die join("\t", ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub ) ) , "\n"; 
	}
	
	if ($m2 == $m1 and $m1 > 1 ){
		$snp_base = "MANNUAL";
		$per = "NA";
		return ( "MANNUAL", "MANNUAL" );
	}
	if( $per == 0 ){
		return ("NONE", "NONE");
	}else{
		return  ("$ref_sub" . "=>" . $snp_base, $per );
	}	
}



sub remove_arrow {#( \$str_modified,  $arrow_num )
	my ( $str_ref,  $num) = @_;
	
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
		if( $string =~ /\^/){
			my $index = $-[0];
			
		#	if($debug){ print STDERR "before:$string\n"}
			substr($string, $index, 2, "\$");
		#	if($debug){ print STDERR "after::$string\n"}
			
		}else{
			die $string;
		}
	}
	
	 if (  $string =~ /\^/){
		print STDERR "moved_string:\t", $string, "\n";
		print STDERR "raw_string:\t", $$str_ref, "\n";
		print $num, "\n\n";
		die;
	}
	$$str_ref = $string;
}

sub remove_dollar{ # ( \$str_modified,  $arrow_num  )
	my ( $str_ref,  $num) = @_;
	
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
		if( $string =~ /\$/){
			my $index = $-[0];
		#	if($debug){ print STDERR "before:$string\n"}
			substr($string, $index, 1, "");
		#	if($debug){ print STDERR "after::$string\n"}
			
		}else{
			die $string;
		}
	}
	
	 if (  $string =~ /\$/){
		print STDERR "moved_string:\t", $string, "\n";
		print STDERR "raw_string:\t", $$str_ref, "\n";
		print STDERR $num, "\n\n";
		die;
	}
	$$str_ref = $string;
}
sub remove_insersion_chars{
	my ($str_ref,  $num, $chr_sub, $pos_sub, $i_sub) = @_;
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
	#	if($debug){ print STDERR "before:$string\n"}
	#	if( $string =~  s/[,\.]\+[0-9]+[ACGTNacgtn]+//){
		if( $string =~  /[,\.]?\+([0-9]+)/){
			my $index = $-[0];
			my $num_of_bp = $1;
	#		if($debug){ print STDERR "after::$string\n"}
		#	if($debug and $pos_sub == 880237){ print STDERR $string, "\t", $index, "\n"}
			 if ( $num_of_bp >= 10){
				#print STDERR "larger than 9\n\n";
				
				die $pos_sub if ($num_of_bp >= 100 );
				
				substr($string, $index,  $num_of_bp + 4, ""); 
			 }else{
				substr($string, $index,  $num_of_bp + 3, ""); 
			 }
		}else{
			print STDERR "$chr_sub; $pos_sub; ", $string, "; $i_sub\n";
			die ;
		}
	}
	
	if ($string =~/\+/){
		print STDERR "moved_string:\t", $string, "\n";
		print STDERR "raw_string:\t", $$str_ref, "\n";
		
		print STDERR "num: $num\n";
		
		die $string, "\n", $$str_ref , "\n\n";

	}
	
	$$str_ref = $string
	
}

sub remove_deletion_chars {
	my ($str_ref,  $num, $chr_sub, $pos_sub, $i_sub) = @_;
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
		#if( $string =~ s/[,\.]?-[0-9]+[ACGTNacgtn]+//){
		#12592042
		if( $string =~  /[,\.]?-([0-9]+)/){
			my $index = $-[0];
			my $num_of_bp = $1;
	#		if($debug){ print STDERR "after::$string\n"}
		#	if($debug and $pos_sub == 12592042){ print STDERR $string, "\t", $index, "\n"}
			 if ( $num_of_bp >= 10){
				#print STDERR "larger than 9\n\n";
				
				die $pos_sub if ($num_of_bp >= 100 );
				substr($string, $index,  $num_of_bp + 4, ""); 
			 }else{
				substr($string, $index,  $num_of_bp + 3, ""); 
			 }
		}else{
			print STDERR "$chr_sub; $pos_sub; ", $string, "; $i_sub\n";
			die ;
		}
	}
	
	if ($string =~/\-/){
		print STDERR "moved_string:\t", $string, "\n";
		print STDERR "raw_string:\t", $$str_ref, "\n";
		
		print STDERR "num: $num\n";
		
		die $string, "\n", $$str_ref , "\n\n";

	}
	
	$$str_ref = $string
}

sub remove_char_and_one_after{
	my ( $str_ref, $char , $num) = @_;
	
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
		if( $string =~ /$char/){
			my $index = $-[0];
			substr($string, $index, 2, "");
			
		}else{
			die $string;
		}
	}
	
	 if (  $string =~ /$char/){
		print STDERR "moved_string:\t", $string, "\n";
		print STDERR "raw_string:\t", $$str_ref, "\n";
	}
	$$str_ref = $string;
}


sub max_sort{
	my @a = @_;
	my @b = sort {$b <=> $a} @a;
	return @b;
}

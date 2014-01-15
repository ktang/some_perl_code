#!/usr/bin/perl -w

use strict;
# /Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/mpileup/513loci_downstream/self_percentage_cutoff_Jun13/db/func_out_uniq_ln/1_WT_0_SNP_only_AA2_func_uniq.txt
#chr	pos	ref	dep	seq	codon_change	AA_change	gene	strand	func
#chr1	2044397	G	20	AAaA,AaAA,AaaAAAAAaa	GTG=>ATG	V=>M	AT1G06670	+	nuclear DEIH-boxhelicase


my $usage = "\n$0 <input> <output>\n\n";
die $usage unless (@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if (-e $output);
open (OUT, ">$output") or die;
open(IN, $input) or die;

my $h = <IN>;
chomp $h;
my @h_a = split "\t", $h;
#print join("\t", (@h_a[0..2], "percentage", @h_a[3..$#h_a]) );

#my $type = $h_a[-1];
#$type =~ s/qual/mut/;

print OUT join("\t", (@h_a[0..2], "type", "percetage", @h_a[3..$#h_a]) ), "\n";

while (<IN>){
	chomp;
	my @a = split "\t";
	my ( $chr, $pos, $ref ) = @a[0..2];
	
	my $dep = $a[3];
	my $seq = $a[4];
	
	my ($type, $per) = categorize($chr, $pos, $ref , $dep, $seq,  1);
	if ($type eq "NO") {
		print  join("\t", ("NO", @a[0..2], $type , @a[3..6])), "\n";
		next;
	}
#	print OUT join("\t", (@a[0..2], $type , @a[3..6])), "\n";
	print OUT join("\t", (@a[0..2], $type, $per, @a[3..$#a])), "\n";
}
close IN;
close OUT;

exit;
sub categorize{
	my ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub) = @_;
	if ($dep == 0){
		#return "NA";
		die ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub) ;
	}
	
	my $str_modified = $str;
	
	my $arrow_num  = ( $str_modified =~ tr/^/^/ );
	if($arrow_num > 0){
		remove_arrow( \$str_modified,  $arrow_num );
	}

	my $plus  = ( $str_modified =~ tr/+/+/ );
	if($plus > 0){
		remove_insersion_chars(\$str_modified,  $plus, $chr_sub, $pos_sub, $i_sub);
	}
	my $minus = ( $str_modified =~ tr/-/-/ );
	if ( $minus > 0 ){
		remove_deletion_chars( \$str_modified, $minus, $chr_sub, $pos_sub, $i_sub);
	}
	
	
	my $dot_num   = ( $str_modified =~ tr/././  );
	my $comma_num = (  $str_modified =~ tr/,/,/   );
	my $A_num = ( $str_modified =~ tr/Aa/Aa/  );
	my $C_num = ( $str_modified =~ tr/Cc/Cc/  );
	my $G_num = ( $str_modified =~ tr/Gg/Gg/  );
	my $T_num = ( $str_modified =~ tr/Tt/Tt/  );

	my $N_num    = ( $str_modified =~ tr/Nn/Nn/  );
	my $star_num = ( $str_modified =~ tr/*/*/  );
	
	my $sum =  $dot_num + $comma_num + $A_num + $C_num + $T_num + $G_num  + $plus + $minus +$N_num + $star_num ;
	
	if($sum != $dep){
		print STDERR join("\t", ($i_sub, $chr_sub, $pos_sub, $ref_sub, $dep, $str,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num  , $plus , $minus ,$N_num , $star_num)), "\n";
	}
	
	my ($m1, $m2, $m3, $m4, $m5, $m6) = max_sort($A_num, $C_num, $G_num, $T_num, $plus, $minus + $star_num );
	
	
	my $max = $m1;
	my $hit_num = 0;
	my $snp_base = "";
	
#	my $per = sprintf( "%.2f" , $max / $dep * 100);
	if( $max == $plus  		 ) {$hit_num ++; $snp_base = "INS"; }
	if( $max == $minus + $star_num   ) {$hit_num ++; $snp_base = "DEL"; }
	
	if( $max == $A_num 		 ) {$hit_num ++; $snp_base  = "A"; }
	if( $max == $C_num 		 ) {$hit_num ++; $snp_base .= "C"; }
	if( $max == $G_num 		 ) {$hit_num ++; $snp_base .= "G"; }
	if( $max == $T_num 		 ) {$hit_num ++; $snp_base .= "T"; }
	
	
	if ($max == 1) {
	#	return "NO";#code
	}
	
	
	unless ($hit_num == 1) {
		print STDERR "$i_sub hit_num:$hit_num\n";
		if ( $max != $minus + $star_num and  $max !=  $plus) {
			print STDERR join("\t", ($i_sub, $snp_base, $chr_sub, $pos_sub, $ref_sub, $dep, $sum, $str, $str_modified  ,
					 $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus , $N_num , $star_num )), "\n";
			#return ("MANUAL");#code
		}
		
	}
	
	my $per = sprintf("%.2f", 100 - ($dot_num+$comma_num)/$dep *100 );
	
	return  ("$ref_sub" . "=>" . $snp_base,$per);
	
	
}



sub remove_arrow {#( \$str_modified,  $arrow_num )
	my ( $str_ref,  $num) = @_;
	
	my $string = $$str_ref;
	
	for  ( my $i = 1;$i <= $num ; $i++ ){
		if( $string =~ /\^/){
			my $index = $-[0];
			
		#	if($debug){ print STDERR "before:$string\n"}
			#substr($string, $index, 2, "\$");
			substr($string, $index, 2, "");
		#	if($debug){ print STDERR "after::$string\n"}
			
		}else{
			#die $string;
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
	#	if( $string =~  /[,\.]?\+([0-9]+)/){
		if( $string =~  /\+([0-9]+)/){
			my $index = $-[0];
			my $num_of_bp = $1;
	#		if($debug){ print STDERR "after::$string\n"}
		#	if($debug and $pos_sub == 880237){ print STDERR $string, "\t", $index, "\n"}
			 if ( $num_of_bp >= 10){
				#print STDERR "larger than 9\n\n";
				
				die $pos_sub if ($num_of_bp >= 100 );
				#G-1A SNP and indel
			#	substr($string, $index,  $num_of_bp + 4, ""); 
				substr($string, $index - 1,  $num_of_bp + 4, ""); 
			 }else{
			#	substr($string, $index,  $num_of_bp + 3, ""); 
				substr($string, $index - 1,  $num_of_bp + 3, ""); 
			 }
		}else{
			print STDERR "$chr_sub; $pos_sub; ", $string, "; $i_sub\n";
		#	die ;
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
	#	if( $string =~  /[,\.]?-([0-9]+)/){
		if( $string =~  /-([0-9]+)/){
			my $index = $-[0];
			my $num_of_bp = $1;
	#		if($debug){ print STDERR "after::$string\n"}
		#	if($debug and $pos_sub == 12592042){ print STDERR $string, "\t", $index, "\n"}
			 if ( $num_of_bp >= 10){
				#print STDERR "larger than 9\n\n";
				
				die $pos_sub if ($num_of_bp >= 100 );
			#	substr($string, $index,  $num_of_bp + 4, ""); 
				substr($string, $index - 1,  $num_of_bp + 4, ""); 
			 }else{
			#	substr($string, $index,  $num_of_bp + 3, ""); 
				substr($string, $index - 1,  $num_of_bp + 3, ""); 
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
	
	$$str_ref = $string;
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

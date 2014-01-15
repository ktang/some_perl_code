#!/usr/bin/perl -w

#~/misc/Zhu_Xiaohong/downstream/dep5/step1_all_common_and_other_23:38:06_N=618$ less not_all_common_3404loci.txt 

#chr1    6324    T       23      T=>INS  100.00  23      13      T=>INS  100.00  13      19      T=>INS  94.74   18      6       T=>INS  66.67   4       11      T=>INS  90.91   10      14      T=>INS  92.86   13      .+1A.+1      1 A,+1a.+1A,+1a.+1A,+1a.+1A.+1A.+1A.+1A.+1A.+1A.+1A,+1a.+1A.+1A.+1A,+1a.+1A,+1a.+1A,+1a    .+1A.+1A.+1A.+1A,+1a,+1a.+1A,+1a,+1a,+1a,+1a.+1A.+1A    a.+1A.+1A.+1A,+1a,+1a.+1A,+1a.+1A.+1A,+1a,+1a,+1a,+1a.+1A,+1a,+1a.+1A,+1a
#				  4 	+4		
			# 4 + 5*4 = 24		
use strict;

my $debug = 0;

my %lables;
my %nums;
#~/misc/Zhu_Xiaohong/downstream/dep5/step1_all_common_and_other_23:57:04_N=622$ less not_all_common_3404loci.txt  | perl ../../src/count_mut_type_number.pl > not_all_common_3404loci_mutation_type_count.txt 

while (<>){
	chomp;
	my @a = split "\t";
#	my $indel = 0;
	my ( $chr, $pos, $ref ) = @a[0..2];
	
	my $x = 0;
	
	for  (my $i = 4; $i<=24; $i+=4){
		$x ++;
		
		my $lab = $a[$i];
		
		$lables{ $lab } = 1;
		
		$nums{$x}->{$lab} ++;
		
		
	}
	
	
}

foreach my $lab ( sort keys %lables){
	print $lab, "\t";
	for my $x (1..5){
		my $n = 0;
		if (defined $nums{$x}->{$lab}  ){
			$n = $nums{$x}->{$lab};
		}
		print $n, "\t";
	}
	
	my $x = 6;
	my $n = 0;
	if (defined $nums{$x}->{$lab}  ){
		$n = $nums{$x}->{$lab};
	}
	print $n, "\n";
}

exit;

sub categorize{
	my ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub) = @_;
	my $str_modified = $str;
	
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
	
	my $N_num    = ( $str =~ tr/Nn/Nn/  );
	my $star_num = ( $str =~ tr/\*/\*/  );
	
	my $sum =  $dot_num + $comma_num + $A_num + $C_num + $T_num + $G_num  + $plus + $minus +$N_num + $star_num ;
	
#	if( $dep !=  $sum and $sum - $dep != 1){ # $arrow_num + $dollar_num
		#print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $str, $str_modified  ,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus,  $arrow_num , $dollar_num  )), "\n", 
	#	 print STDERR join("\t", ( $chr_sub, $pos_sub, $ref_sub, $dep, $sum, $str, $str_modified  ,  $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num, $plus , $minus , $N_num , $star_num )), "\n", 
#	}
	
	my ($m1, $m2, $m3, $m4, $m5, $m6) = max_sort($A_num, $C_num, $G_num, $T_num, $plus, $minus + $star_num );
	
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
		return ( "MANNUAL", "NA", 0);
	}
	if( $per == 0 ){
		return ("NONE", "NONE", 0);
	}else{
		return  ("$ref_sub" . "=>" . $snp_base, $per, $max );
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

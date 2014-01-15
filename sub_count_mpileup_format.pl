#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

close(IN);
close(OUT);

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}

sub categorize{
	my ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub) = @_;
#	my ($str) = @_;
	if ($dep == 0){
		return "NA";
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
	
#	my $N_num    = ( $str =~ tr/Nn/Nn/  );
#	my $star_num = ( $str =~ tr/*/*/  );

	my $N_num    = ( $str_modified =~ tr/Nn/Nn/  );
	my $star_num = ( $str_modified =~ tr/*/*/  );
	
	my $sum =  $dot_num + $comma_num + $A_num + $C_num + $T_num + $G_num  + $plus + $minus +$N_num + $star_num ;
	
	if($sum != $dep){
		print STDERR join("\t", ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub, $dot_num , $comma_num , $A_num , $C_num , $T_num , $G_num  , $plus , $minus ,$N_num , $star_num)), "\n";
	}
	
	
	if($dot_num + $comma_num == $dep){
		return "REF";
	}
	
	my ($m1, $m2, $m3, $m4, $m5, $m6) = max_sort($A_num, $C_num, $G_num, $T_num, $plus, $minus + $star_num );
	
	
	my $max = $m1;
	my $snp_base = "NA";
	

	if   ( $max == $A_num 		 ) { $snp_base = "A"; }
	elsif( $max == $C_num 		 ) { $snp_base = "C"; }
	elsif( $max == $G_num 		 ) { $snp_base = "G"; }
	elsif( $max == $T_num 		 ) { $snp_base = "T"; }
	elsif( $max == $plus  		 ) { $snp_base = "INS"; }
	elsif( $max == $minus + $star_num) { $snp_base = "DEL"; }
	else{
		#die join("\t", ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub ) ) , "\n"; 
		print STDERR join("\t", ($chr_sub, $pos_sub, $ref_sub, $dep, $str, $i_sub ) ) , "\n"; 
	}
	
	return ( $snp_base);

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


sub max_sort{
	my @a = @_;
	my @b = sort {$b <=> $a} @a;
	return @b;
}
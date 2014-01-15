#!/usr/bin/perl -w

use strict;
# 0	 1	2	3		4						5	 6		7						8
#chr2    990     N       15      ^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A   444444444444444 17      AAAA^'A^'A^'A^'A^'A^'A^'A^-A^'A^'A^'A^'A^'A     33334444444444444

#my $usage = "\n$0 <input> <col_number_of_mpileup_seq_1based> STDOUT\n\n";
#die $usage unless (@ARGV == 2);
my $usage = "\n$0 <input> STDOUT\n\n";
die $usage unless (@ARGV == 1);

my $input = shift or die;
#my $col_num = shift or die;
my $col_num = 5;

die unless (-e $input);
open(IN, $input) or die;

#chr	pos	ref	dep	seq	codon_change	AA_change	gene	strand	func
#chr1	9260067	G	19	,$.,.,aaA,AAaa,.A,.A	GCT=>ACT	A=>T	AT1G26770	+	expansin A10


my $h = <IN>;
my @h_a = split "\t", $h;
print join("\t", (@h_a[0..2], "percentage", @h_a[3..$#h_a]) );

my %r;

while (<IN>){
	my @a = split "\t";
	#print if ($a[0] eq "chr");
	my ( $chr, $pos, $ref ) = @a[0..2];
	my $dep = $a[$col_num - 2];
	my $seq = $a[$col_num - 1];
	my $per = categorize($chr, $pos, $ref , $dep, $seq, 1);
	#print print join("\t", (@a[0..2], $per, @a[3..$#a] )) ;
	$r{$chr}->{$pos} = join("\t", (@a[0..2], $per, @a[3..$#a] )) ;
	
}
close IN;

foreach my $chr (sort keys %r) {
	foreach my $pos (sort {$a<=>$b} keys %{$r{$chr}}){
		print $r{$chr}->{$pos};
	}
}

exit;


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
	
	my $per_sub = sprintf("%.2f", 100* ($dot_num + $comma_num) /$dep );
	
	return (100-$per_sub);
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

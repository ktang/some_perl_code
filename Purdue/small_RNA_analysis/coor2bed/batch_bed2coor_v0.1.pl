#!/usr/bin/perl -w

#v0.1
#output comment before head

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <indir> <outdir>\n\n";
die $usage unless(@ARGV == 2);

my $indir = shift or die;
my $outdir = shift or die;

die unless (-d $indir);
die unless (-d $outdir);

opendir(DIR, $indir) or die;
my @files = grep /\.txt$/, readdir DIR;
closedir(DIR);

foreach my $file(@files){
	if($file =~ /(\S+)\.txt$/){
		my $pre = $1;
		my $input = File::Spec->catfile($indir, $file);
		my $output = File::Spec->catfile($outdir, $pre . "_coor.txt");

		die unless (-e $input);
		die if( -e $output);

		open(IN, $input) or die "cannot open $input: $!";

		die if(-e $output);
		open(OUT, ">$output") or die "cannot open $output: $!";
		
		my $h;
		while($h = <IN>){
			if($h =~ /^#/){
				print OUT $h;
			}else{
				last;
			}
		}
		
	#	my $h = <IN>;
		chomp ($h);
		my @h_a = split "\t", $h;
		print OUT join("\t", ( "coor", @h_a[3..$#h_a]) ) ,"\n";
		#print OUT join("\t", ("chr", "start", "end") ), "\n";

		#my ($chr, $s, $e) ;

		while (<IN>) {
			chomp;
			my @a=split "\t";
			next if ($a[0] =~ /[PMt]/ or $a[0] eq "chrC");
			my $chr = simple_chr($a[0]);

			print OUT join ("\t", ("$chr:$a[1]-$a[2]", @a[3..$#a]) ),"\n";
	
		}


		close(IN);
		close(OUT);
	}
	else{
		die $file;
	}
	
}

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
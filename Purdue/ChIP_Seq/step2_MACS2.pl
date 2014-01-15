#!/usr/bin/perl -w

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $genome_size = "1.2e8";
my $keep_dup = "all";

my $usage = "\n$0 \n\n <treatment_ChIP> <control_Input> <outdir> <outpre> \n\n";
die $usage unless(@ARGV == 4);

my $tfile     = shift or die;
my $cfile     = shift or die;
my $outdir    = shift or die;
my $outpre    = shift or die;

my $log_file = File::Spec->catfile($outdir, $outpre . "_MACS2_log.txt");
die "log file exist\n\n" if(-e $log_file);

my $cmd = "time macs2 callpeak -t $tfile -c $cfile -g $genome_size --keep-dup $keep_dup --outdir $outdir -n $outpre 2>> $log_file";

print STDERR $cmd, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd, "\n";
	close OUT;
	`$cmd`;
}
exit;

=head
if(!$debug){
	die unless(-e $);
	die if(-e $);
}

print STDERR $cmd_, "\n\n";
unless($debug){
	open(OUT, ">>$log_file") or die;
	print OUT $cmd_, "\n\n";
	close OUT;
	`$cmd_`;
}



=cut

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
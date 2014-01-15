#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;



#print Kai_Module::round( 10.8), "\n\n";
#print Kai_Module::simple_chr("chr4"), "\n\n";

print round( 10.8), "\n\n";
print simple_chr("chr4"), "\n\n";


my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}

=head
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
=cut

exit;

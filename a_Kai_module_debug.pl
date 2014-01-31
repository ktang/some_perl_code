#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);
# my ($volume,$directories,$file) =          File::Spec->splitpath( $path );


BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

use utf8;#可以吗？
use strict;
use File::Spec;

my $x = "5.88,5.88,11.76,312";

my @a = split ",", $x;

print Kai_Module::get_ind_of_maximum(\@a), "\n";

#print Kai_Module::round( 10.8), "\n\n";
#print Kai_Module::simple_chr("chr4"), "\n\n";

#print round( 10.8), "\n\n";
#print simple_chr("chr4"), "\n\n";
#Kai_Module::first_call("ashgdk");
# 14 ,,...-2TA.-4TATA*.-2TA.-2TA,-2ta.-2TA.-2TA,-2ta.-2TA HCI=DA!DD?DD@D
=head
my $usage = "$0 \n <dep> <pileup> <qual>\n\n";
die $usage unless (@ARGV == 3);

my $dep = shift or die;
my $pileup  = shift or die;
my $qual= shift or die;
my $usage = "$0 \n <input>\n\n";
my $input = shift or die;

open(IN, $input) or die;
while (<IN>) {
    chomp;
    my ($dep,  $pileup , $qual) = split "\t";
    my @a = Kai_Module::split_pileup($dep,  $pileup , $qual);
    print STDERR join("\t", @a), "\n";

}

=cut

#my @a = Kai_Module::split_pileup(14, ",,...-2TA^\$.-4TATA*.-2TA.-2TA,-2ta.-2TA.-2TA,-2ta.-2TA", "HCI=DA!DD?DD\@D");
#print STDERR join("\t", @a), "\n";

#my @a = Kai_Module::split_pileup($dep,  $pileup , $qual);
#print STDERR join("\t", @a), "\n";


exit;

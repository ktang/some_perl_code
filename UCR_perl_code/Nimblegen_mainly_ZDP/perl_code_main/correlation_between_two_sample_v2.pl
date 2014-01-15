#!/usr/bin/perl -w


=head
This program mainly task is to input two files of the two bilogical replicates,
finally get the correlation coefficient and the mean of difference
=cut

##begin of main

#input

my $usage = "$0 <input_file_1> <input_file_2>";

die $usage unless (@ARGV == 2);

$ARGV[0] =~ /(.+)\.gff/;
my $pre1 = $1;
$ARGV[1] =~ /(.+)\.gff/;
my $pre2 = $1;
my $out = join('_',$pre1,$pre2);

#print "$output\n";


system("
cut -f 6 $ARGV[0] > temp_1.txt;
cut -f 6 $ARGV[1] > temp_2.txt;

paste -d '\t' temp_1.txt temp_2.txt > R_in;
rm temp_1.txt;
rm temp_2.txt;
echo 'R...';
R < ./correlation_v2.R --no-save;
rm R_in;
#mv correlation.txt ${$ARGV[1].$ARGV[2]}.cor;
#mv mean_of_diff.txt ${$ARGV[1].$ARGV[2]}.mean;
paste name.txt out.txt > $out.stat;
");

print "\a\n\n";

exit;
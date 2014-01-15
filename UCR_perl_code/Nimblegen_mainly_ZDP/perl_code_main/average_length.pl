#!/usr/bin/perl


=head
This program mainly task is to input TAMALg results  to
finally get the everage length of the peaks.
=cut

##begin of main

#input

my $usage = "$0 <input_file_1>";

die $usage unless (@ARGV == 1);

$ARGV[0] =~ /(\S+)L2/;
$pre = $1;

print "\n\n\n\n$pre\n\n\n";
$out = join("",$pre,"avar_len_of_peaks");

system("
cut -f 4,5 $ARGV[0] > R_in;
echo \"R...\";
R --no-save < ./average_length.R;
rm -f R_in;
mv R_out $out;
");

print "\n\n\n$out\n\n";
#!/usr/bin/perl


=head
This program mainly task is to input two files of the two bilogical replicates,
finally get the correlation coefficient and the mean of difference
=cut

##begin of main

#input

my $usage = "$0 <input_file_1> <input_file_2>";

die $usage unless (@ARGV == 2);

`
cut -f 6 $ARGV[0] > temp_1.txt;
cut -f 6 $ARGV[1] > temp_2.txt;
paste -d '\t' temp_1.txt temp_2.txt > R_in;
echo "R...";
R --no-save < ./correlation.R;
`;
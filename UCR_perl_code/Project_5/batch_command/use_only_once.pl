#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 17 Nov 2010

=cut

use strict;

my $adapter1 = "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";

my $adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

my $in1 = "/Users/kaitang/Desktop/zhu_lab_data_orignal/100724_HWUSI_EAS1501_0006_61RRRAXX";
my $out1 = "/Users/kaitang/Desktop/Project_5/soap_output/100724_HWUSI_EAS1501_0006_61RRRAXX";


my $in2 = "/Users/kaitang/Desktop/zhu_lab_data_orignal/100822_HWUSI-EAS1501_00007_61RF3AAXX_PE1";
my $out2 = "/Users/kaitang/Desktop/Project_5/soap_output/100822_HWUSI-EAS1501_00007_61RF3AAXX_PE1";

my $cmd_157_1 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in1/157_1PairedEnd $out1/157_1PairedEnd $adapter1 $adapter2 C24";

my $cmd_30_19 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in1/30_19PairedEnd $out1/30_19PairedEnd $adapter1 $adapter2 C24";

my $cmd_330_15 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in1/330_15Ros1PairedEnd $out1/330_15Ros1PairedEnd $adapter1 $adapter2 C24";

my $cmd_Mut77 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in1/Mut77PairedEnd $out1/Mut77PairedEnd $adapter1 $adapter2 Col0";

my $cmd_col_seq5 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in2/Col_gl1/fastq/seq5 $out2/Col_gl1/seq5 $adapter1 $adapter2 C24";

my $cmd_col_seq6 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in2/Col_gl1/fastq/seq6 $out2/Col_gl1/seq6 $adapter1 $adapter2 C24";

my $cmd_ros2_seq7 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in2/ros2/fastq/seq7 $out2/ros2/seq7 $adapter1 $adapter2 C24";

my $cmd_ros2_seq8 = "time /Users/kaitang/Desktop/perl_code/Project_5/batch_command/PE_fastq_batch_SOAP_pileup.pl $in2/ros2/fastq/seq8 $out2/ros2/seq8 $adapter1 $adapter2 C24";

print "$cmd_157_1\n";
`$cmd_157_1`;

print "$cmd_30_19\n";
`$cmd_30_19`;

print "$cmd_330_15\n";
`$cmd_330_15`;

print "$cmd_Mut77\n";
`$cmd_Mut77`;

print "$cmd_col_seq5\n";
`$cmd_col_seq5`;


print "$cmd_col_seq6\n";
`$cmd_col_seq6`;

print "$cmd_ros2_seq7\n";
`$cmd_ros2_seq7`;

print "$cmd_ros2_seq8\n";
`$cmd_ros2_seq8`;

exit;
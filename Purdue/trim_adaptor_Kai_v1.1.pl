#!/usr/bin/perl -w
# Copyright (C)  2010-
# Program:			trim_adaptor_Kai.pl
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		Jan 23, 2011
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	Feb 4, 2011
# Description:		modified from Dr. Renyi Liu's (UCR) program
#**************************
# Version: 1.0 input fastq, output fasta
# Version: 1.1 input fastq, output fastq
#**************************
# e-mail:tangkai.ustc@gmail.com

my $version="1.1";
use strict;
use warnings;

use Getopt::Long;
use FindBin;
use Bio::SeqIO;

my $debug = 0;
my $i = 0;
### Command line parameters -------------------------------------------------------------------
my $input	= "";	#input fastq file;
my $output	= "";	#output fasta file;
my $adaptor_3	= "";	#adaptor string;
my $minLen		= "";
my %CMD;
GetCom();

print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";



my $adaptor_1st_1nt = substr($adaptor_3, 0, 1);
my $adaptor_1st_2nt = substr($adaptor_3, 0, 2);
my $adaptor_1st_3nt = substr($adaptor_3, 0, 3);
my $adaptor_1st_4nt = substr($adaptor_3, 0, 4);
my $adaptor_1st_5nt = substr($adaptor_3, 0, 5);
my $adaptor_1st_6nt = substr($adaptor_3, 0, 6);
my $adaptor_1st_7nt = substr($adaptor_3, 0, 7);
my $adaptor_1st_8nt = substr($adaptor_3, 0, 8);
my $adaptor_1st_9nt = substr($adaptor_3, 0, 9);
my $adaptor_1st_10nt = substr($adaptor_3, 0, 10);
my $adaptor_1st_11nt = substr($adaptor_3, 0, 11);
my $adaptor_1st_12nt = substr($adaptor_3, 0, 12);
my $adaptor_1st_13nt = substr($adaptor_3, 0, 13);

my $adaptor_1st_13nt_rev = reverse $adaptor_1st_13nt;

my ($trim_1,$trim_2, $trim_3, $trim_4, $trim_5, $trim_6, $trim_7,
	$trim_8, $trim_9, $trim_10, $trim_11, $trim_12, $trim_13_over, $no_trim) =
	(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); # number of seqs subject to different trim
my $no_short_seqs = 0;	

print STDERR "Trim 3' adaptor for file $input\n";

my $seqin	= Bio::SeqIO->new(-file=>$input, -format=>'fastq');
#my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fasta',width=>200);
############################# 1.1
my $seqout	= Bio::SeqIO->new(-file=>">$output", -format=>'fastq',width=>200);
##############################


while(my $seq = $seqin->next_seq){
		my $seqstr = $seq->seq;
		my $new_seqstr = $seqstr;
		my $flag_no	= 0;	my $flag_ge13	= 0;
		
		
		if($seqstr =~ /(\w*$adaptor_1st_13nt)\w*/){
			
			$trim_13_over++;
			$new_seqstr = $1;
			$flag_ge13	= 1;
			$new_seqstr = reverse $new_seqstr;
			if ($new_seqstr =~ /\w*($adaptor_1st_13nt_rev\w*)/){
				$new_seqstr = $1;
			}
			$new_seqstr = reverse $new_seqstr;
#			if ($new_seqstr =~ /$adaptor_1st_13nt/){
				$new_seqstr =~ s/$adaptor_1st_13nt//;
#			}
			
			if($debug){
				$i ++;
				if ($i < 2){
					print STDERR "Doing seq ", $seq->id, " ", $seq->seq, "\n";
					print STDERR "match >= 9nt; new_seqstr: ", $new_seqstr, "\n";
					my	$debug_new_seqstr = reverse $new_seqstr;
					print STDERR "rev new_seqstr: ", $debug_new_seqstr, "\n";
					$debug_new_seqstr  =~ s/$adaptor_1st_13nt_rev/lol/;
					print STDERR "adaptor_rev replaced: ", $debug_new_seqstr, "\n";
				}
			}

		}elsif(substr($seqstr, length($seqstr) - 12) eq $adaptor_1st_12nt){
			$trim_12++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 12);
			
		}elsif(substr($seqstr, length($seqstr) - 11) eq $adaptor_1st_11nt){
			$trim_11++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 11);
			
		}elsif(substr($seqstr, length($seqstr) - 10) eq $adaptor_1st_10nt){
			$trim_10++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 10);
			
		}elsif(substr($seqstr, length($seqstr) - 9) eq $adaptor_1st_9nt){
			$trim_9++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 9);
			
		}elsif(substr($seqstr, length($seqstr) - 8) eq $adaptor_1st_8nt){
			$trim_8++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 8);
			
		}elsif(substr($seqstr, length($seqstr) - 7) eq $adaptor_1st_7nt){
			$trim_7++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 7);
			
		}elsif(substr($seqstr, length($seqstr) - 6) eq $adaptor_1st_6nt){
			$trim_6++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 6);
			
		}elsif(substr($seqstr, length($seqstr) - 5) eq $adaptor_1st_5nt){
			$trim_5++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 5);
			
		}elsif(substr($seqstr, length($seqstr) - 4) eq $adaptor_1st_4nt){
			$trim_4++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 4);
			

		}elsif(substr($seqstr, length($seqstr) - 3) eq $adaptor_1st_3nt){
			$trim_3++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 3);
			
		}elsif(substr($seqstr, length($seqstr) - 2) eq $adaptor_1st_2nt){
			$trim_2++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 2);
			
		}elsif(substr($seqstr, length($seqstr) - 1) eq $adaptor_1st_1nt){
			$trim_1++;
			$new_seqstr = substr($seqstr, 0, length($seqstr) - 1);
			
		}else{
			$no_trim++;	$flag_no	= 1;
			
		}

#		if($new_seqstr =~ /N/){
#			$no_seqs_has_N++;
#			next;
#		}

		if(length($new_seqstr) < $minLen){
			$no_short_seqs++;
		}
		else{
			
###########################1.1			
			
		
#			my $index = length($new_seqstr);
#			my $sub_qual = $seq->subqual(1,$index);
#			print STDERR "index:$index\n";
#			my @a = ((@{$qual})[0..$index]);
#			my @b = @a[0..$index];
		#	my @a_b = @{$qual}[0..$index];
#			$seq->qual(\@a);
			
	#		print STDERR scalar(@$qual),"\n";
#			print STDERR scalar(@a_b),"\n";
#			print STDERR "new_qual:",scalar(@a),"\nb:",scalar(@b),"\n";
		#	print STDERR join(' ',@new_qual),"\n";
#			my @sub_qual = @{seq->subqual(1,length($new_seqstr))};
#			$seq->qual(\@sub_qual);
=head2 qual
The length of the returned value always matches the length
           of the sequence.
 Title   : qual
 Usage   : $qual_values  = $obj->qual($values_string);
 Function:

           Get and set method for the meta data starting from residue
           position one. Since it is dependent on the length of the
           sequence, it needs to be manipulated after the sequence.

           The length of the returned value always matches the length
           of the sequence.

 Returns : reference to an array of meta data
 Args    : new value, string or array ref or Bio::Seq::PrimaryQual, optional

=cut
###########################1.1				
			$seq->seq($new_seqstr);
###########################1.1					
			my $sub_qual = $seq->subqual(1,length($new_seqstr));
			$seq->qual($sub_qual);
###########################1.1					
			$seqout->write_fastq($seq);
		}
		
}		
	

summary();
sub_end_program();
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}

############################################################################################################
######################                  summary
############################################################################################################
sub summary{
	print join("\t",(">=13","12","11","10","9", "8", "7", "6", "5", "4", "3", "2", "1", "0","short than $minLen bp")), "\n";
	print join("\t",($trim_13_over,$trim_12,$trim_11,$trim_10,$trim_9, $trim_8, $trim_7, $trim_6, $trim_5, $trim_4, $trim_3, $trim_2, $trim_1, $no_trim, $no_short_seqs)),"\n";
}

############################################################################################################
######################                  Read command line parameters
############################################################################################################
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory:
--adaptor	STRING
--input		STRING	fastq_file
--output	STRING	fastq_file
--minLen	number	minium length of seq can be left after triming
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "adaptor=s", "input=s", "output=s","minLen=n");
	
	# Mandatory params
	die("Please specify adaptor sequence\n") unless defined($CMD{adaptor});
	die("Please specify input file\n") unless defined($CMD{input});
	die("Please specify output file\n") unless defined($CMD{output});
	die("Please specify minLen\n")	unless defined($CMD{minLen});
	$adaptor_3	= $CMD{adaptor};	#adaptor string;
	$input	= $CMD{input};	#input fastq file;
	$output= $CMD{output};	#output fasta file;
	$minLen		= $CMD{minLen};
	print STDERR "adaptor:$adaptor_3\ninput:$input\noutput:$output\nminium length is $minLen\n";
}


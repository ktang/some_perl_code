#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

# get command-line arguments, or die with a usage statement
# my $usage = "$0 \n <infile> <infileformat> <outfile> <outfileformat> \n\n";
#die $usage unless (@ARGV == 4);

#    /Users/tang58/DataBase/TAIR10/Proteins/TAIR10_pep_20101214

my $usage = "$0 \n <infile> \n\n";
die $usage unless (@ARGV == 1);
my $infile = shift or die $usage;
#my $infileformat = shift or die $usage;

#create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new( '-file' => "<$infile",
                              '-format' => "fasta");

# write each entry in the input file to the output file
#  while (my $inseq = $seq_in->next_seq) {
#    $seq_out->write_seq($inseq);
# }

#/Users/tang58/DataBase/TAIR10/Proteins/TAIR10_pep_20101214
my $inseq = $seq_in->next_seq;

print join("\t", ("id", $inseq-> id)), "end\n";
print join("\t", ("accession_number", $inseq-> accession_number)), "end\n";
print join("\t", ("desc", $inseq->desc )), "end\n";
print join("\t", ("primary_id", $inseq->primary_id )), "end\n";
print join("\t", ("alphabet", $inseq-> alphabet)), "end\n";
print join("\t", ("object_id", $inseq-> object_id)), "end\n";
print join("\t", ("namespace", $inseq->namespace )), "end\n";
print join("\t", ("display_name", $inseq-> display_name)), "end\n";
print join("\t", ("description", $inseq->description )), "end\n";
print join("\t", ("annotation", $inseq->annotation )), "end\n";
#print join("\t", ("", $inseq-> )), "end\n";
#print join("\t", ("", $inseq-> )), "end\n";


exit;
#id	AT1G51370.2end
#accession_number	unknownend
#desc	| Symbols:  | F-box/RNI-like/FBD-like domains-containing protein | chr1:19045615-19046748 FORWARD LENGTH=346end
#primary_id	AT1G51370.2end
#alphabet	proteinend
#object_id	unknownend
#namespace	end
#display_name	AT1G51370.2end
#description	| Symbols:  | F-box/RNI-like/FBD-like domains-containing protein | chr1:19045615-19046748 FORWARD LENGTH=346end
#annotation	Bio::Annotation::Collection=HASH(0x100a48ac8)end

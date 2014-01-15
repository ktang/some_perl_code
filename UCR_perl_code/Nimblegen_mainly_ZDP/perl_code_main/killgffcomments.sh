#!/bin/sh


#killgffcomments.sh
#takes a gff file.
#FILE MUST END WITH .gff
#takes one file on command line - outputs automatically to a new file with "nopound" in the name


awk '{if ($1 !~ /\#/) {print $0;};}' $1 > $2


#!/bin/sh

#sortgfffilev03.sh
#produces a sorted gff file
#INPUT
#takes one inputfile on command line and the output filename on the command line
#EXAMPLE
#	 sortgfffilev03.sh 98498.gff first.gff

awk  '{OFS="\t";if ($1 !~ /\#/) {print $1 OFS (sprintf ("%0.10u", $2)) OFS (sprintf ("%0.10u", $3))};};' $1 > myzfhold.gff
sort -k1,2 myzfhold.gff > ${2}
rm myzfhold.gff
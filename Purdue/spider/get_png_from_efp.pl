#!/usr/bin/perl -w
#right:
#http://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&show_obsolete=F&name=at4g18650&sub_type=gene&SEARCH_EXACT=4&SEARCH_CONTAINS=1

#wrong:
#http://www.arabidopsis.org/servlets/Search?type=general&search_action=&method=&show_obsolete=&name=&sub_type=gene&SEARCH_EXACT=4&SEARCH_CONTAINS=1&words=AT4G18650
use strict;
use WWW::Mechanize;
use Getopt::Long;

my $debug = 0;
#efp-hRdXMM.png
#efp-4xZ_jw.png

# http://bbc.botany.utoronto.ca/efp/cgi-bin/output/efp-QAwtuT.png
# http://bbc.botany.utoronto.ca/efp/cgi-bin/output/efp-hRdXMM.png
# http://bbc.botany.utoronto.ca/efp/cgi-bin/output/efp-QAwtuT.png
#http://affymetrix.arabidopsis.info/narrays/search.pl?f1=1&s1=ATGE_21
#http://bbc.botany.utoronto.ca/efp/cgi-bin/output/efp-GkYMVv.png
my $search = shift or die "Must specify a search term";
my $w = WWW::Mechanize->new;
# http://bbc.botany.utoronto.ca/efp/cgi-bin/efpWeb.cgi?primaryGene=AT4G18650&modeInput=Absolute

# $w->get( "http://www.arabidopsis.org/" );
 $w->get( "http://bbc.botany.utoronto.ca/efp/cgi-bin/efpWeb.cgi?primaryGene=" . $search. "&modeInput=Absolute" );
$w->success or die "Can't read the search page!\n";

if ($w->content =~ m/loadPopup.+src="output\/efp-(\S+)\.png.*/){
			#my $imgurl = "http://bbc.botany.utoronto.ca/efp/cgi-bin/output/$1.png";
			my $imgurl = "http://bar.utoronto.ca/efp_arabidopsis/cgi-bin/output/efp-" . $1 . ".png";
			if($debug) {
				print STDERR $imgurl, "\n\n";
				print STDERR $1, "\n";
			}
				
			my $filename = "$search.png";
# Bring down the image if we don't already have it
			if ( -e $filename ) { warn "Already have $filename\n"; }
			else {

				warn "Saving $filename...\n";
				$w->get( $imgurl, ":content_file"=>$filename );

			}
}
exit;

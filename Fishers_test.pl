#!/usr/bin/perl -w

#         word2   ~word2
#  word1    n11      n12 | n1p
# ~word1    n21      n22 | n2p
#          --------------
#           np1      np2   npp

#/Volumes/Macintosh_HD_2/idm1_new_met_data/script/fishers_exact_test
use strict;

#use Text::NSP::Measures::2D::Fisher2::left;
use Text::NSP::Measures::2D::Fisher2::twotailed;
#use Text::NSP::Measures::2D::Fisher::twotailed;
my $npp = 37; my $n1p = 27; my $np1 = 20;  my $n11 = 14;

my  $left_value = calculateStatistic( n11=>$n11,
                                      n1p=>$n1p,
                                      np1=>$np1,
                                      npp=>$npp);

my $errorCode;

if( ($errorCode = getErrorCode()))
  {
    print STDERR $errorCode." - ".getErrorMessage();
  }
  else
  {
    print getStatisticName." value(Fisher) for bigram is ".$left_value, "\n";
  }

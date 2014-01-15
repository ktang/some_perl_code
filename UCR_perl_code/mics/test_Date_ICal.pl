#! /usr/bin/perl -w
use strict;
use Test::Simple tests => 2;

use Date::ICal;

my $ical = Date::ICal -> new;         # create an object
ok( defined $ical );                # check that we got something
ok( $ical->isa('Date::ICal') );     # and it's the right class

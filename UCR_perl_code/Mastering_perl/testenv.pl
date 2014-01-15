#!/usr/bin/perl

print "Content-type: text/plain\n\n" if $ENV{REQUEST_METHOD};

foreach my $key (sort keys %ENV)
{
	printf "%-20s %s\n", $key,$ENV{$key};	
}
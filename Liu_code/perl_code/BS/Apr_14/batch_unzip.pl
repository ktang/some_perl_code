#!/usr/bin/perl -w
use strict;

opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		print STDERR "unzipping under $dir\n";
		chdir($dir);
		opendir(SDIR, ".");
		my @files = grep {/\.tar\.gz/} readdir SDIR;
		
		foreach my $file(@files){

		    `tar zxf $file`;
					}
	    chdir('..');
	}
}

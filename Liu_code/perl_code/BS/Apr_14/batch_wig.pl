#!/usr/bin/perl -w
use strict;

opendir(CDIR, ".");
my @dirs = readdir CDIR;
foreach my $dir(@dirs){
    if((-d $dir) && ($dir =~ /lane/)){
		print STDERR "creating wig file under $dir\n";
		chdir($dir);
		opendir(SDIR, ".");
		my @files = grep {/_forw\.txt/} readdir SDIR;
		#my $out = "conversion_error.txt";
		foreach my $file(@files){
		    my $pre = "NONE";
		    if($file =~ /(s_\d+)/){
			    $pre = $1;
		    }
			my $rev_file = $pre . "_rev.txt";
			
			my $out = $pre . "_methylation.wig";
		    `perl ../create_wig_file.pl $pre $file $rev_file > $out`;
			`gzip $out`;
		}
	    chdir('..');
	    
	 
	}
}

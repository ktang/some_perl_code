#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
my $debug = 0;
my $usage = "$0 <do>";
die $usage unless(@ARGV == 1 and $ARGV[0] eq "do");

#my $cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021424/SRR052365/SRR052365.sra";
#print STDERR $cmd, "\n\n";

if(!$debug){
#	`$cmd`;
}

# ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP002%2FSRP002100/SRR037710/SRR037710.sra

my %pres = ( 
					
				"SRR037791" => "SRR037791.sra",
				"SRR037792" => "SRR037792.sra",
				"SRR037879" => "SRR037879.sra",
				"SRR037880" => "SRR037880.sra",
				"SRR037793" => "SRR037793.sra",
				
				"SRR037794" => "SRR037794.sra",
				"SRR037710" => "SRR037710.sra",
				"SRR037788" => "SRR037788.sra",
				"SRR037785" => "SRR037785.sra",
				"SRR037786" => "SRR037786.sra",
				
				"SRR066816" => "SRR066816.sra",
				"SRR066817" => "SRR066817.sra"
			);

foreach my $pre( sort keys %pres){
	
	my $file = $pres{$pre};
	die "$file" if(-e $file);
	
	die unless ($file =~ /$pre/);
	
	
	my $cmd = "time wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP002%2FSRP002100/" . $pre . "/" . $file;
	print STDERR $cmd, "\n\n";

	if(!$debug){
		`$cmd`;
	}
}

exit;
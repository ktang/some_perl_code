#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
my $debug = 0;
my $usage = "$0 <do>";
die $usage unless(@ARGV == 1 and $ARGV[0] eq "do");

my $cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021424/SRR052365/SRR052365.sra";
print STDERR $cmd, "\n\n";

if(!$debug){
	`$cmd`;
}

#ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021424/SRR052370/SRR052370.sra
foreach my $i(67..70){
	#my $cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021424/SRR052365/SRR052365.sra";
	my $cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021424/SRR0523" . $i . "/SRR0523" .$i .".sra";
	print STDERR $cmd, "\n\n";

	if(!$debug){
		`$cmd`;
	}
}

exit;
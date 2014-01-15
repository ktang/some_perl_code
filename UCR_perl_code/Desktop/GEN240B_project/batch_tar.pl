#!/usr/bin/perl -w
# batch download human chromosome from
# ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/


for (my $i = 2; $i <= 25; $i++)
{
	if ($i <= 22)
	{
		my $cmd = "gunzip Homo_sapiens.GRCh37.57.dna.chromosome.$i.fa.gz";
		`$cmd`;	
	}	
	elsif($i == 23)
	{
	  	my $cmd = "gunzip Homo_sapiens.GRCh37.57.dna.chromosome.X.fa.gz";
		`$cmd`;	
	}
	elsif($i == 24)
	{
	  	my $cmd = "gunzip Homo_sapiens.GRCh37.57.dna.chromosome.Y.fa.gz";
		`$cmd`;	
	}
	elsif($i == 25)
	{
	  	my $cmd = "gunzip Homo_sapiens.GRCh37.57.dna.chromosome.MT.fa.gz";
		`$cmd`;	
	}
	
}

  print "\a";
  exit;

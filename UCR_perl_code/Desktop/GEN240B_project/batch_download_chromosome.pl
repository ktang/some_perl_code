#!/usr/bin/perl -w
# batch download human chromosome from
# ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/


for (my $i = 3; $i <= 25; $i++)
{
	if ($i <= 22)
	{
		my $cmd = "wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.57.dna.chromosome.$i.fa.gz";
		`$cmd`;	
	}	
	elsif($i == 23)
	{
	  	my $cmd = "wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.57.dna.chromosome.X.fa.gz";
		`$cmd`;	
	}
	elsif($i == 24)
	{
	  	my $cmd = "wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.57.dna.chromosome.Y.fa.gz";
		`$cmd`;	
	}
	elsif($i == 25)
	{
	  	my $cmd = "wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.57.dna.chromosome.MT.fa.gz";
		`$cmd`;	
	}
	
}

  print "\a";
  exit;

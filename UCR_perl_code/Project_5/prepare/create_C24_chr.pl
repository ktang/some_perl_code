#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 17 Nov 2010

only work for 5 chrs.
not ChrM,ChrC

Nov 17,2010
=cut

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $debug=0;

my $usage = "$0 <Col-0 Chr> <C24_SNP> <output>";
die $usage unless(@ARGV == 3);

my ($Col, $C24_SNP,$output) = @ARGV[0..2];

my $chr_in = Bio::SeqIO -> new ('-file' => "<$Col",
			    '-format' => "fasta");
			    
my $out = Bio::SeqIO->new(-file=>">$output", -format=>"fasta");

			    
open (SNP, "$C24_SNP")
	or die "cannot open SNP file:$!";
	
my @SNPs = <SNP>;

close (SNP);

my %chrs;

while (my $chr = $chr_in->next_seq)
{
	my $chr_name = $chr->id;
	
	my $chr_number = -1;
	if ($chr_name =~ /Chr(\d)/)
	{
		$chr_number = $1;
		if($debug)
		{
			print "chr_number=$chr_number\n"
		}
	}
	if ($chr_number != -1)
	{
		$chrs{$chr_number} = $chr->seq();	
		if ($debug)
		{
			print STDERR "chromosome:$chr_number\t$chrs{$chr_number}\n";	
		}
	}
}

$chr_in->close;

for(my $i = 0; $i <= $#SNPs; $i++)
{
	my $this = $SNPs[$i];
	
	if($debug)
	{
		print STDERR "in for i:$i\n$this";	
	}
	
	chomp $this;
	
	my @pts = split "\t",$this;
	my ($chro,$position,$ref_base,$sub_base)= @pts[1..4];
	if ($debug)
	{
			print STDERR "$this\n$chro,$position,$ref_base,$sub_base\n";
	}
	
	my $orignal = substr($chrs{$chro},$position-1,1);
	if ($orignal eq $ref_base)
	{
		substr($chrs{$chro},$position-1,1,$sub_base);	
	}
	else 
	{
		print STDERR "warning:",$this,"\n";
		print STDERR "in fasta is $orignal\n";
#		die;	
	}
}
  
  foreach my $num (sort {$a<=>$b} keys %chrs)
  {
  	
  	if ($debug)
  	{
  		print STDERR "num=$num\n";	
  	}
  	my $seq_name = "Chr".$num;
  	my $seq = Bio::Seq->new(-seq => $chrs{$num}, -id=>$seq_name);
  	$out->write_seq($seq);
  }
  
   $out->close;
  
 
exit;

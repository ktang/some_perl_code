#!/usr/bin/perl -w

use strict;

#START HERE MAIN PROGRAM

my $usage = "$0 <Ecker_data> <peaks_gff>";
die $usage unless (@ARGV == 2);

$ARGV[0] =~ /(\S+)\.gff/;
my $pre_e = $1;
my $out_e = join ("",$pre_e,"_sorted.gff");

$ARGV[1] =~  /(\S+)\.gff/;
my $pre_p = $1;
my $out_p =  join ("",$pre_p,"_sorted.gff"); 

system ("sortgfffilev03.sh $ARGV[0] sorted.temp");
system("killgffcomments.sh sorted.temp $out_e");
system ("sortgfffilev03.sh $ARGV[1] sorted.temp");
system("killgffcomments.sh sorted.temp $out_p");
system("rm -f sorted.temp");
my $ecker_file = $out_e;
my $peakfilenm=$out_p;
my %nums;



#load gff lines into @lines
open(ECKER,"<$ecker_file");
my @eckers=<ECKER>;
close(ECKER);

#load peak file
open(PEAKF,"<$peakfilenm");
my @peaks=<PEAKF>;
close(PEAKF);

my $c_num = $#eckers + 1; #note that this is also the sequence length
my $p_num = $#peaks + 1;
my $total_overlap = 0;

# **CALCULATE VALUE FOR EACH PEAK
my $slopval=0;
my ($p,$c) = (0 , 0);

LOOP: while ( ($p < $p_num) and ($c < $c_num))
{
	my $thispeak = $peaks[$p];
	my $thisC = $eckers[$c];
	chomp $thispeak;
	chomp $thisC;
	my @pts_p = split "\t", $thispeak;
	my @pts_c = split "\t", $thisC;
	if( $pts_p[0] eq $pts_c[0])
	{
		if( ($pts_c[3] >= $pts_p[3] )and($pts_c[3] <= $pts_p[4] ) )
		{
			$nums{$pts_p[0]}->{$pts_p[3]}++;
			$c++;
			$total_overlap++;
			next LOOP;		
		}
		elsif ($pts_c[3] < $pts_p[3] ) 
		{
			$c++;
			next LOOP;
		}
		elsif($pts_c[3] > $pts_p[4] )
		{
			$p++;
			next LOOP;
		}
	}

	elsif ( $pts_p[0] lt $pts_c[0] ) 
	{
		$p++;
		next LOOP;
	}
	
	elsif ( $pts_p[0] gt $pts_c[0] )
	{
		$c++;
		next LOOP;
	}

}

foreach my $chr (sort keys %nums)
{
	foreach my $pos (sort {$a <=> $b} keys %{$nums{$chr}})
	{
		print "$chr\t$pos\t$nums{$chr}->{$pos}\n";
	}
}

print STDERR  "\a";
exit;

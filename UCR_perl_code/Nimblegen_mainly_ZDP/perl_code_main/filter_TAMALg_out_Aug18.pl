#!/usr/bin/perl -w

=head2 filter_TAMALg_Aug18.pl

PROGRAM: filter_TAMALg_Aug18.pl
GOAL: filter the output of TAMALg to make sure that the value 
      of mutant is positive before substraction.
INPUT: <TAMALg_output.gff> <mutant_normalized_value_gff>
file must be gff format

OUTPUT:
TAMALg_output_filtered.gff


=cut


sub isoverlap
{

	my $acor=$_[0];
	my $bcor=$_[1];
	my $ccor=$_[2];
	my $dcor=$_[3];
	
	if ($acor>=$ccor && $acor<=$dcor)
	{return 1;}
	if ($bcor>=$ccor && $bcor<=$dcor)
	{return 1;}
	if ($ccor>=$acor && $ccor<=$bcor)
	{return 1;}
	if ($dcor>=$acor && $dcor<=$bcor)
	{return 1;}
	return 0;
}


my $usage = "$0 <TAMALg_output.gff> <mutant_normalized_value_gff>";

die $usage unless (@ARGV == 2);

my $fname1=$ARGV[0]; #TAMALg filenam
my $fname2=$ARGV[1]; #mutant mean filename


my ($mend,$mstart,$pend,$pstart);
my ($m,$p) = (0,0);

$fname1 =~ /(\S+)\.gff/;

my $pre = $1;
my $out = $pre."_filter.gff";
my $sorted = $pre."_sorted.gff";

open (OUT ,">$out")
 or die("cannot open output file $out:$!");


system ("sort -k1,4 $fname1 > $sorted");

open (SHORT,"<$sorted");
my @peak_file=<SHORT>;
close(SHORT);



open (LONG,"<$fname2");
my @mean_file = <LONG>;
close(LONG);

 


my @flag;

for(my $i = 0; $i <= $#peak_file; $i++ )
{
	$flag[$i] = 1;	
}


 
while ( ( $m <= $#mean_file ) and ( $p <= $#peak_file ) )
{
	my $l_m = $mean_file[$m];
	my $l_p = $peak_file[$p];
	my @pts_m = split "\t", $l_m;
	my @pts_p = split "\t", $l_p;
	$mstart = $pts_m[3];
	$mend = $pts_m[4];
	$pstart = $pts_p[3];
	$pend = $pts_p[4];
	if ( ( $pts_m[0] eq $pts_p[0] ) and  ( isoverlap($mend,$mstart,$pend,$pstart) ) )
	{
		if ( $pts_m[5] < 0 )
		{
			$flag[$p] = -1;
			print "$mean_file[$m], $peak_file[$p]";
			$p++;	
		}
	}
	
	elsif ( ( $pts_m[0] eq $pts_p[0] ) and ( $mend < $pstart  ) )
	{
		$m++;	
	}
	
	elsif ( ( $pts_m[0] eq $pts_p[0] ) and ( $mstart > $pend  ) )
	{
		$p++;	
	}
	
	else 
	{
		if( $pts_m[0] lt $pts_p[0] )
		{
			$m++;	
		}
		elsif( $pts_m[0] gt $pts_p[0] )
		{
			$p++;
		}
	}
}



for (my $i = 0; $i <=  $#peak_file; $i++ )
{
	if( $flag[$i] == 1)
	{
		print OUT $peak_file[$i];
	}	
}

close (OUT);
print "\a";
exit;

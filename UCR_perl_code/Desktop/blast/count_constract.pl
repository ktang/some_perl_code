#!/usr/bin/perl -w
# count output of blast
use strict;

my $usage = "$0 <in>";
die $usage unless(@ARGV == 1);

open (IN,"$ARGV[0]")
 or die "cannot open:$!";

my $total = 0;
my $line;
my @parts;
my %genes;

while ($line = <IN>)
{
	chomp $line;
	@parts = split /\t/, $line;
	if ($parts[0] eq $parts[1] && !(exists $genes{$parts[0]}))
	{
		$genes{$parts[0]} = $parts[2];
		$total++;
	}
	
}

print $total,"\n";

close(IN);

$total = 0;

open (IN,"$ARGV[0]")
 or die "cannot open:$!";
 while ($line = <IN>)
{
	chomp $line;
	@parts = split /\t/, $line;
	if ($parts[0] ne $parts[1] )
	{
		if( !(exists $genes{$parts[0]})	)
		{
			$total++;
			print $line,"\n";
		}
		elsif ($parts[2] >= $genes{$parts[0]})
		{
			$total++;	
			delete $genes{$parts[0]};
			print $line,"\n";
		}
	}
	
}

close(IN);

print "contrast:",$total,"\n";
$total = 0;
foreach (keys %genes)
{$total++;}
print $total,"\n";
exit;
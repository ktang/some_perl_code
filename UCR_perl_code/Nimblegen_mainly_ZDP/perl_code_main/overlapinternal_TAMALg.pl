#!/usr/bin/perl

=head2 overlapinternal.pl

PROGRAM: overlapinternal.pl
GOAL: looks for overlaps in a single file; outputs the file with merged pieces
INPUT: one filename then a slopvalue (set to zero for strict overlaps)
file must be gff format

OUTPUT:
to STDOUT

EXAMPLE:
overlapinternal.pl  E2F1_MARKMAX_H_T05P05_S50_CHR_G2.gff 0

SPEEDNOTE:
this could be sped up A LOT by sorting the files to begin, and then stopping the looping as soon
as the first field doesn't match the field



=cut


sub isoverlappingslop
{

	my $extrabit=$_[0]; #allows for sloppiness - if this is set to zero, then there is no sloppiness
	my $acor=$_[1]-$extrabit;
	my $bcor=$_[2]+$extrabit;
	my $ccor=$_[3]-$extrabit;
	my $dcor=$_[4]+$extrabit;
	
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


$fname1=$ARGV[0]; #input filenam
#$fname2=$ARGV[1]; #second filename
$slopval=$ARGV[1] + 0; #extra bit to add to each piece to allow near overlaps

open (DOG,"<$fname1");
@first=<DOG>;
close(DOG);


=head2 older

open (DOG2,"<$fname2");
@second=<DOG2>;
close(DOG2);




$first_numlines=$#first +1;
$second_numlines=$#second+1;



#print "first numlines = $first_numlines second numlines = $second_numlines\n";

for ($i=0; $i<$first_numlines; $i++)
{
	$first_overlaps[$i]=0;
}

for ($i=0; $i<$second_numlines; $i++)
{
	$second_overlaps[$i]=0;
}

=cut

=head2 junk

for ($i=0; $i < $first_numlines; $i++)
{
	print "i is $i firstval is $first[$i]\n";
}

=cut

#note use $i as index here
$i=0;

while ($#first>0)
{
	$thisone=shift(@first);
	@pieces = split "\t", $thisone;
	$fstart=$pieces[3];
	$fend=$pieces[4];
	$overlapfound=0; #flag for overlap
	LOOKLOOP: for ($j=0; $j <= $#first; $j++)
	{
		@spieces = split "\t", $first[$j];
		$sstart=$spieces[3];
		$send=$spieces[4];
		#print "fstart is $fstart fend is $fend sstart is $sstart send is $send\n";
		
		
		if (($pieces[0] eq $spieces[0]) and isoverlappingslop($slopval, $fstart,$fend,$sstart,$send))
		{
			$overlapfound=1;

			#GET NEW PARAMS FOR LINE
			if ($fstart<$sstart)
				{$tstart=$fstart;}
			else
				{$tstart=$sstart;}
			if ($fend>$send)
				{$tend=$fend;}
			else
				{$tend=$send;}
			@newpiece=@pieces;
			$newpiece[3]=$tstart;
			$newpiece[4]=$tend;
			$mergeline=join "\t", @newpiece;

			#REFORM ARRAY (DELETE ELEMENT, UNSHIFT MERGED ONE
			@tempcp=(@first[0..($j-1)],@first[($j+1)..($#first)]);
			unshift  @tempcp, $mergeline;
			@first=@tempcp;

			#END LOOP AND RESTART. Note that we readded the piece, so we'll just work with it more

			last LOOKLOOP;

			#$second_overlaps[$j] = $second_overlaps[$j] + 1;
			#print $first[$i] . $second[$j] . "******\n";
		}
		#print "second overlaps 0 is $second_overlaps[0]\n";
		#print "i is $i j is $j\n";
		
	}
if ($overlapfound==0)
	{print $thisone;}

#reset $overlapfound
$overlapfound=0;
}

#PRINT THE LAST ONE!
print $first[0];


=head2 old code

#print "second overlaps[9] is " .  $second_overlaps[9] . "\n";
$fcount=0;
$scount=0;

for ($i=0; $i < $first_numlines; $i++)
{
	if ($first_overlaps[$i]>0)
	{
		$fcount++;
	}
}

$fper=sprintf "%2.2f",(100*$fcount/$first_numlines);

for ($i=0; $i < $second_numlines; $i++)
{
	if ($second_overlaps[$i]>0)
	{
		$scount++;
	}
}

$sper=sprintf "%2.2f",(100*$scount/$second_numlines);

#OUTPUT HERE
#added in overlappercent_oneline.pl
#format: slopval, firstfilename, firstoverlapped number, total items in firstfile, percent overlapped, secondfilename, numoverlapped, totalitems in second, percent

print "$slopval\t$fname1\t$fcount\t$first_numlines\t$fper\t$fname2\t$scount\t$second_numlines\t$sper\n";


=cut


=head2 older output format

print "$fname1 vs $fname2 with S of $slopval\n";
print "$fname1 $fcount / $first_numlines  for $fper %\n";
print "$fname2 $scount / $second_numlines  for $sper %\n";
#print "$fname2 length is $second_numlines num_overlapped: $scount percent: $sper\n";

=cut


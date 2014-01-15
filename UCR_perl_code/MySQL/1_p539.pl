#!/usr/bin/perl -w

use strict;
use DBI;

my $dbh = DBI->connect("DBI:mysql:host=localhost;database=mysql","root");

my $sth = $dbh->prepare("select * from lucy.wage");
$sth->execute();
while (my $ref = $sth->fetchrow_hashref())
{
	print ($ref->{monthly});
}
$sth->finish();

$dbh->disconnect();

use strict;
use warnings;

use Gene1;

print "Object 1:\n\n";
my $obj1 = Gene1->new(
	name		=> "Aging",
	organism	=> "Homo sapiens",
	chromosome	=> "23",
	pdbref		=> "pdb999.net"
);

print $obj1 ->name ,"\n";
print $obj1 ->organism ,"\n";
print $obj1 ->chromosome ,"\n";
print $obj1 ->pdbref ,"\n";

print "Object 2:\n\n";
my $obj2 = Gene1->new(
	organism	=> "Homo sapiens",
	name		=> "Aging",
);

print $obj2 ->name ,"\n";
print $obj2 ->organism ,"\n";
print $obj2 ->chromosome ,"\n";
print $obj2 ->pdbref ,"\n";

print "Object 3:\n\n";
my $obj3 = Gene1->new(
	organism	=> "Homo sapiens",
	chromosome	=> "23",
	pdbref		=> "pdb999.net"
);

print $obj3 ->name ,"\n";
print $obj3 ->organism ,"\n";
print $obj3 ->chromosome ,"\n";
print $obj3 ->pdbref ,"\n";

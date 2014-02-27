package NAME;


use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.0.0;
@ISA         = qw(Exporter);
@EXPORT      = qw( round);
@EXPORT_OK   = qw( round);
%EXPORT_TAGS = (
		DEFAULT => [qw(&round)]
	       );
#%EXPORT_TAGS = ( DEFAULT => [qw(&func1)],
 #                Both    => [qw(&func1 &func2)]);

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}


1;
__END__

BEGIN { push @INC, '/Users/tang58/scripts_all/perl_code/Modules' }
use Kai_Module;

if ($debug) {
	print STDERR , "\n";
}

http://www.perlmonks.org/?node_id=102347

First we get a namespace by declaring a package name. This helps ensure our module's
functions and variables remain separate from any script that uses it.

Use strict is a very good idea for modules to restrict the use of global variables.
See use strict warnings and diagnostics or die for more details.

We need to use the Exporter module to export our functions from the MyModule::
namespace into the main:: namespace to make them available to scripts that 'use' MyModule.

We pacify strict with the use vars declaration of some variables. We can use an 'our' declaration in 5.6+

We now set a $VERSION number and make Exporter part of MyModule using the @ISA. See perlboot for all the
gory details on what @ISA is or just use it as shown.

@EXPORT contains a list of functions that we export by default, in this case nothing. Generally the less
you export by default using @EXPORT the better. This avoids accidentally clashing with functions defined
in the script using the module. If a script wants a function let it ask.

@EXPORT_OK contains a list of functions that we export on demand so we export &func1 &func2 only if
specifically requested to. Use this in preference to just blindly exporting functions via @EXPORT.
You can also export variables like $CONFIG provided they are globals not lexicals scoped with my
(read declare them with our or use vars).

%EXPORT_TAGS. For convenience we define two sets of export tags. The ':DEFAULT' tag exports only &func1;
the ':Both' tag exports both &func1 &func2. This hash stores labels pointing to array references.
In this case the arrays are anonymous.

We need the 1; at the end because when a module loads Perl checks to see that the module returns
a true value to ensure it loaded OK.
You could put any true value at the end (see Code::Police) but 1 is the convention.


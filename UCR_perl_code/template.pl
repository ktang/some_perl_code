#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 

Nov 17,2010
=cut

use strict;
my $usage = "$0 <input_bam>";
die $usage unless(@ARGV == 1);

exit;


#!/usr/bin/perl -w
# Copyright (C)  2010-
# Program:			
# Author:			Kai Tang <tangkai.ustc@gmail.com>
# Program Date:		
# Modifier:			Kai Tang <tangkai.ustc@gmail.com>
# Last Modified:	
# Description:	
#**************************
# Version: 1.0
#**************************
# e-mail:tangkai.ustc@gmail.com

my $version="1.0";
print STDERR ("\n==================| $0 start |==========================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); #ÔËÐÐ¿ªÊ¼Ê±¼ä
print STDERR "Now = $Time_Start\n\n";


############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #Ê±¼ä×Ó³ÌÐò
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	my $end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}
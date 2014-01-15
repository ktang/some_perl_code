#!/usr/bin/perl -w

# check whether probe is ordered
# by the position in chr

use strict;

my $usage = "$0 <input_gff>";

die $usage unless(@ARGV == 1);

open (IN,  "<$ARGV[0]")
  or die "Cannot open input_File: $!";
  
  
  my @lines;
  
  @lines = <IN>;
  close(IN);

 
  my @ps1;
  my $last_chr = "chr";
  my ($last_pos);
  for(my $i = 0; $i <= $#lines; $i++)
  {
  		@ps1 = split "\t", $lines[$i];
  		
  		if ($last_chr eq "chr"){
  			$last_chr = "chr1";
  			$last_pos = $ps1[3] + 0;	
  		}
  		elsif ($ps1[0] eq $last_chr){
  			my $this_pos = $ps1[3] + 0;
  			if ($this_pos < $last_pos){
  				print "$lines[$i-1]$lines[$i]";
  				die; 
  			}
  			else {
  				$last_pos = $ps1[3] + 0;
  			}
  		}
  		else {
  			$last_chr = $ps1[0];
  			$last_pos = $ps1[3] + 0;
  		}
  }
  
  print STDERR "\a";
  exit;
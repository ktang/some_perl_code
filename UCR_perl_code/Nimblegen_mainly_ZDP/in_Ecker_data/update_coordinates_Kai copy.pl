#!/usr/bin/perl
use warnings;
use strict;

my $file = shift;

if (!defined($file) || $file eq "" ){
    print "To run the script try: \n\n";
    print "./update_coordinates.pl <assembly_changes.txt> <file_to_update.gff> <new_file_name (opt)>\n";
    exit;
}
my $gff = shift;
if (!defined($gff) || $gff eq "" ){
    print "To run the script try: \n\n";
    print "./update_coordinates.pl <assembly_changes.txt> <file_to_update.gff> <new_file_name (opt)>\n";
    exit;
}

my $outputfile = shift;
if (!defined($outputfile) || $outputfile eq ""){
    print "WARNING: output file will be named tair9_" . $gff . "\n";
    $outputfile = "tair9_" . $gff;
}

main();

sub main {
    open (OUT, ">$outputfile");
    for ( my $i=1; $i<8; $i++ ){
        my %insert_hash;
        my %delete_hash;
        
        #
        # read in the file with all the insertions and deletions
        #
        open(FILE, $file);
        while (my $line = <FILE>){
            chomp($line);
            my @lines = split("\t", $line);
            if ($lines[0] == $i){
                if ($lines[2] eq "insertion" && @lines == 4) {
                    $insert_hash{ $lines[1] } = length($lines[3]);
                } elsif ($lines[2] eq "insertion" && @lines == 5){
                    $insert_hash{ $lines[1] } = $lines[3];
                } elsif ($lines[2] eq "deletion" && @lines == 4){
                    $delete_hash{ $lines[1] } = length($lines[3]);
                } elsif ($lines[2] eq "deletion" && @lines == 5){
                    $delete_hash{ $lines[1] } = $lines[3];
                }
            }
        }
        close(FILE);
    
        #
        # read in the other file and start the transform
        # 
        open (FILE, $gff);
        while (my $line = <FILE>){
            chomp($line);
            my @lines = split("\t", $line);
            my $chr = "CHR" . $i;
            if ($i == 6) { $chr = "CHRC"; } elsif ($i==7) { $chr = "CHRM"; }
            my $start = $lines[3];
            my $end = $lines[4];
            my $cur_chr = uc($lines[0]);
	    my $warn = 0;
            if ( uc($cur_chr) eq $chr || uc($cur_chr) eq $i ) {
                foreach my $key ( keys %delete_hash ){
                    
                    # check deletion falls within the region
                    if ( ($lines[3] >= $key) && ($lines[4] < ($key + $delete_hash{$key}) ) ){
                        $start = ($key - 1);
                        $end = $key;
                        print $lines[0]."\t".$lines[1]."\t".$lines[2]."\t".$start."\t".$end."\t".$lines[5]."\t".$lines[6]."\t".$lines[7]."\t".$lines[8]."\t";
                        print "#WARNING:ObsoleteFeature;TAIR8_old_coordinates:$chr:$key:" . ($delete_hash{$key}+$key) . ":+\n";
                       
#mask Kai
# print OUT "#WARNING:ObsoleteFeature;TAIR8_old_coordinates:$chr:$key:" . ($delete_hash{$key}+$key) . ":+\n";
                        $warn = 1;
		    # check to see if 5' end falls in a deletion
                    } elsif ( ($lines[3] >= $key ) && ($lines[3] < ($key + $delete_hash{$key}) ) && ($delete_hash{$key} >= 5) ) {
                        $start = $lines[3];
                    # check to see if 3' end falls in a deletion
                    } elsif ( ($lines[4] >= $key) && ($lines[4] < ($key + $delete_hash{$key}) ) && ($delete_hash{$key} >= 5) ){
                        $end = $lines[4];
                    } else {
                        if ($lines[3] >= $key) { $start = $start - $delete_hash{$key}; }
                        if ($lines[4] >= $key) { $end = $end - $delete_hash{$key}; }
                    }
                }
                foreach my $key ( keys %insert_hash ){
                    if ($lines[3] >= $key) {
                        $start = $start + $insert_hash{$key};
                    }
                    if ($lines[4] >= $key){
                        $end = $end + $insert_hash{$key};
                    }
                }
                print OUT $lines[0]."\t".$lines[1]."\t".$lines[2]."\t".$start."\t".$end."\t";
          #      if ($lines[5]) { print OUT $lines[5]; } print OUT "\t";
           #     if ($lines[6]) { print OUT $lines[6]; } print OUT "\t";
            #    if ($lines[7]) { print OUT $lines[7]; } print OUT "\t";
             #   if ($lines[8]) { print OUT $lines[8]; } 
	#	if ($warn == 1) { print OUT ";WARNING=Deleted Feature"; } print OUT "\t\n";
                print OUT $lines[5]."\t".$lines[6]."\t".$lines[7]."\t".$lines[8]."\n";
            }
        }
        close(FILE);
    }
    close(OUT);
}

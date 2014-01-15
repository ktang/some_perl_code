#!/usr/bin/perl -w

#system("date");
use strict;
use File::Spec;

my $debug = 0;
if($debug){
	print STDERR "debug: $debug\n";
}

# my $script = "/Users/tang58/Weiqiang_idm1_1_Nature_paper/add_annotation_Kai_v1.7_TAIR10.pl";
#my $script = "/Users/tang58/Kai_BS/for_publish/add_annotation_Kai_TAIR10_v1.8.pl";

my $script = "/Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/add_annotation_Kai_TAIR10_v1.8_coordinate.pl";
die unless ( -e $script ); 

my $usage = "\n<indir> <outdir> <pattern>\n\n";

die $usage unless (@ARGV == 3);

my ($indir, $outdir, $pattern) = @ARGV;

die "wrong dir" unless (-d $indir and -d $outdir);

opendir (DIR, $indir) or die "cannot open $outdir: $!";

my @files = grep /$pattern/, readdir DIR;

if ($debug){
	print STDERR join("\n", @files),"\n";
}

foreach my $input(@files){
	if( $input =~ /(\S+)\.txt$/){
		
		my $full_input = File::Spec->catfile($indir, $input);
				
		my $output = File::Spec->catfile( $outdir,  $1 . "_annotated_TAIR10.txt" );
		#die "wrong input" unless (-e "$indir/$input");
		#die "wrong output" unless ( ! (-e "$outdir/$output"));
		#my $cmd = "perl $script $indir/$input > $outdir/$output";
		
		die "wrong input: $full_input" unless (-e $full_input);
		die "wrong output" if  ( -e $output );
		my $cmd = "perl $script $full_input > $output";
		
	
		if($debug){
			print STDERR $cmd, "\n\n";
		}
		else{
			print STDERR "handling $input...";
			`$cmd`;
			print STDERR "\tdone\n";
		}
	}
	else{
		die $input;
	}
}

exit;

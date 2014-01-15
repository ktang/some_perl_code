#!/usr/bin/perl -w

# v0.0
#
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}


#my $usage = "$0 \n <postfix> <outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
#die $usage unless(@ARGV >= 8);

my $usage = "$0 \n <number_of_files> <file> <label> [<file> <label> ...] \n\n";
die $usage unless(@ARGV >= 3);


#my $postfix = shift or die;
#my $outdir = shift or die;
#my $dep_cutoff = shift or die;
my $num_of_files = shift or die;
my $last_index = $num_of_files - 1;

#die unless (-d $outdir);

die unless (@ARGV == 2 * $num_of_files);

my @labels;
my @files;
#my @outputs;

for my $i(0..$last_index){
	$files[$i]  = shift or die;
	$labels[$i] = shift or die;
#	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_isMeth_depth" . $dep_cutoff . "_" . $postfix . ".txt" );
	die $files[$i] unless (-e $files[$i]);
#	die if (-e $outputs[$i]);
}

#print STDERR "depth_cutoff is $dep_cutoff \n\n";
#if($debug){
	print STDERR join("\n", @labels),  "\n\n";
	print STDERR join("\n", @files),   "\n\n";
#	print STDERR join("\n", @outputs), "\n\n";
#}

if($debug){
	exit;
}

my %r_label;
my %r_mut;
#chr     pos     ref     mut_WT_0        per_WT_0        dep_WT_0        seq_WT_0        qual_WT_0       type    type_code
#chr1    2044397 G       G=>A    90.00   20      AAaA,AaAA,AaaAAAAAaa    ACFEJDHICJIIHIJIJEC#    CDS     1
#0	 1	 2	 3	 4	 5
for my $i(0..$last_index){
	my $label = $labels[$i] ;
	open( IN, "<" , $files[$i]) or die "$files[$i]";
	my $h = <IN>;
	while(<IN>){
		chomp;
		my @a = split "\t";
		my ($chr, $pos) = @a[0..1];
		my $mut = $a[3];
		push @{$r_label{$chr}->{$pos}}, $label;
		if(defined $r_mut{$chr}->{$pos}){
			if($r_mut{$chr}->{$pos} ne $mut ){
				print STDERR $r_mut{$chr}->{$pos} , "\t", $_, "\n";
			}
		}else{
			 $r_mut{$chr}->{$pos} = $mut;
		}
	}
	close IN;
}

my %num_h;
foreach my $chr (sort keys %r_label ){
	foreach my $pos (sort {$a<=>$b} keys %{$r_label{$chr}}){
		my $tmp = join("\t", @{$r_label{$chr}->{$pos}});
		$num_h{$tmp}++;
	}
}

foreach my $key (sort keys %num_h){
	print join("\t", ($key, $num_h{$key})),"\n";
}


exit;

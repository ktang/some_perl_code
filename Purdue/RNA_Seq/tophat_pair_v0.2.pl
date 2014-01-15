#!/usr/bin/perl -w

#time tophat2 -p 6
#-G /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf
#-o Fengqiu_WT_thout
# --b2-very-sensitive
#--phred64-quals
#-r 200

#/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome
#CleanData/WT_L1_1.fq.gz CleanData/WT_L1_2.fq.gz 
use strict;
use File::Spec;

my $debug = 1;

print STDERR "incomplete\n\n";
exit;
#my %Phred_h = ( 33=> "--solexa-quals", 64=>"--phred64-quals");
#33 --solexa-quals??
my %Phred_h = ( 33=> " ", 64=>"--phred64-quals");


my $p = 6;
my $r = 200;
my $mode_flag = "--b2-very-sensitive";


my $gtf = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf";
die unless (-e $gtf);

my $genome = "/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome";


#$cmd = "perl $script $sub_indir  $sub_outdir $phred";

#my $usage = "$0 \n<indir> <in1> <in2> <outdir> <phred_33_64>\n\n";
#die $usage unless(@ARGV == 5);

my $usage = "$0 \n<indir>  <outdir> <phred_33_64>\n\n";
die $usage unless(@ARGV == 3);



my $indir = shift or die "indir";
die "indir" unless (-d $indir);
#my $f1 = shift or die "f1";
#my $f2 = shift or die "f2";

my $outdir = shift or die;
die "outdir" if(-d $outdir);

my $phred = shift or die "Phred";
die "Phred = 33 or 64"  unless (defined $Phred_h{ $phred });
my $phred_label = $Phred_h{ $phred };

#my $fq1 = File::Spec->catfile($indir, $f1);
#my $fq2 = File::Spec->catfile($indir, $f2);

#time tophat2 -p 6
#-G /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2012-03-08-18-36-47/Genes/genes.gtf
#-o Fengqiu_WT_thout
# --b2-very-sensitive
#--phred64-quals
#-r 200

#/Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/genome
#CleanData/WT_L1_1.fq.gz CleanData/WT_L1_2.fq.gz

opendir(DIR, $indir ) or die;
my @mates1_name = grep /_R1_\S+\.gz$/, readdir DIR;
closedir DIR;

if($debug){
	print STDERR join("\n", @mates1_name), "\n\n\n";
}
opendir(DIR, $indir ) or die;
my @mates2_name = grep /_R2_\S+\.gz$/, readdir DIR;
closedir DIR;

if($debug){
	print STDERR join("\n", @mates2_name), "\n\n\n";
}

my @mate1_file = add_dir ($indir, \@mates1_name);
my @mate2_file = add_dir ($indir, \@mates2_name);

if($debug){
	print STDERR join("\n", @mate1_file), "\n\n\n";
	print STDERR join("\n", @mate2_file), "\n\n\n";
}

my $in1 = join(",", @mate1_file);
my $in2 = join(",", @mate2_file);

#my $cmd = "tophat2 -p $p -r $r -G $gtf -o $outdir $mode_flag $phred_label $genome $fq1 $fq2";
my $cmd = "tophat2 -p $p -r $r -G $gtf -o $outdir $mode_flag $phred_label $genome $in1 $in2";
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
exit;

sub add_dir{
	my ($dir, $ref) = @_;
	my @a;
	foreach my $f(@{$ref}){
		my $t = File::Spec->catfile($dir, $f);
		push @a, $t;
	}
	return @a;
}

#!/usr/bin/perl -w
# batch map unique seqs to genome
use strict;
use Getopt::Long;

my($indir,$outdir,$ref,$n,$r,$v,$M, $format);

my %CMD;
GetCom();

#my $usage = "$0 <unique seq dir> <soap out dir> <database>";
#die $usage unless(@ARGV >= 3);
#my ($indir, $outdir, $db) = @ARGV[0..2];

my @files = "";
my $ref_file = "";
if ($ref eq "Col0"){ $ref_file = "/Users/tang58/DataBase/TAIR_Col0_genome/index/SOAP/5Chr/5Chr_only_TAIR9_Col0.fas.index";}
elsif($ref eq "C24") {$ref_file = "/Users/tang58/DataBase/C24/index/SOAP/C24_TAIR9_5Chr.fas.index";}
else{die "ref should be C24|Col0"}

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

@files = grep {/\.$format$/} readdir INDIR;

=head
#if ($format eq "fastq"){
#	@files = grep {/\.fastq$/} readdir INDIR;
#}elsif ($format eq "fasta"){
#	@files = grep {/\.fasta$/} readdir INDIR;
#}else{
#	die "file format should be fastq or fasta";
#}
=cut

foreach my $file(@files){
     if($file =~ /(\S+)\.$format$/){
		 my $pre = $1;
		 my $output = $pre . "_vs_"."$ref"."_n$n"."r".$r."v".$v."M".$M.".soapout";
		 my $cmd = "soap -a $indir/$file -D $ref_file -v $v -M $M -r $r -n $n -o $outdir/$output";
		 print $cmd, "\n";
		 #$cmd .= " 2>&1";
		 #my $msg = `$cmd 2>&1`;
		 #print $msg, "\n";
		 `$cmd`;
	 }
}

####################################################################################

sub GetCom{

		my @usage = ("$0
		
	Mandatory:
	--indir		STRING		input_dir_path
	--outdir	STRING		output_dir_path
	--db		STRING		Col0 or C24
	--format	STRING		fasta or fastq
	--n		INT		Filter low quality reads contain more INT bp Ns
	--r		INT		How  to  report repeat hits, 0=none; 1=random one; 2=all
	--v   		INT 		Totally allowed mismatches in one read;at most 2
	--M   		INT 		Match mode for each read or the seed part of read,  which
 						    shouldn't contain more than 2 mismaches
 						      0: exact match only
						      1: 1 mismatch match only
						      2: 2 mismatch match only
						      3: [gap] (coming soon)
						      4: find the best hits	
	\n");
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD,"indir=s","outdir=s","db=s", "format=s","n=i","r=i","v=i","M=i");
	die("Please specify indir\n") unless (defined $CMD{indir});
	die("Please specify outdir\n") unless (defined $CMD{outdir});
	die("Please specify ref\n") unless (defined $CMD{db});
	die("Please specify n\n") unless (defined $CMD{n});
	die("Please specify r\n") unless (defined $CMD{r});
	die("Please specify v\n") unless (defined $CMD{v});
	die("Please specify M\n") unless (defined $CMD{M});
	die("Please specify format\n") unless (defined $CMD{format});
	$indir	=  $CMD{indir};
	$outdir =  $CMD{outdir};
	$ref	=  $CMD{db};
	$n		=  $CMD{n};
	$r		=  $CMD{r};
	$v 		=  $CMD{v};
	$M		=  $CMD{M};
	$format =  $CMD{format};
}
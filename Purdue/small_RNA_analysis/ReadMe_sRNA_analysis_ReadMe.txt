############
#	mapping process
########################
just use TAIR10 genome

################
#	step 1 from unsorted bam to split to 2 bam files  structural RNA (tRNA, rRNA, snRNA, snoRNA) and add NH:i: flag
####################
~/DataBase/TAIR10/GFF_11:38:36_N=511$ less TAIR10_GFF3_genes.gff | grep -E 'tRNA|snRNA|snoRNA|rRNA' | grep Note | less

less TAIR10_GFF3_structural_RNA.txt | perl -F"\t" -lane ' $F[0]=~s/Chr//; if($F[0] eq "M"){$F[0] = "Mt"} if($F[0] eq "C"){$F[0] = "Pt"}  print join("\t", @F)' > TAIR10_GFF3_structural_RNA_simpleChr.txt 

NH:i:1
Number of reported alignments that contains the query in the current record

1 record 

################
#	step2: generate table of 500bp-bin, HNA value (Lee 2012)
################################
#reads_for_normalization:XXX
chr:s-e	hits_number	reads_num	HNA

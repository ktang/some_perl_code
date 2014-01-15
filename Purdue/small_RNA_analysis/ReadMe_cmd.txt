/Volumes/My_Book/20120226_YanJun_sRNA/clean_fa_16:40:42_N=512$ time zcat ../fq_gz/002831_col-1_ATCACG_L001_R1_001.fastq.gz ../fq_gz/002831_col-1_ATCACG_L002_R1_001.fastq.gz |  awk ' { if(NR % 4 == 1 ) { sub( /@/, ">", $1 ) ; print $1 } if(NR % 4 ==2) {print} } ' |  fastx_clipper -a TGGAATTCTCGGGTGCCA -l 18 -c -v -o YanJun_002831_col-1_sRNA_clean_min18.fa


/Volumes/My_Book_2/raw_data_backup/Shanghai/55_RDM16_Col0_background_sRNA/55/result_primary_13:47:03_N=515$ time cat 55.fq | grep TCGTATGCCGTCTTCTGCTTGA | wc -l
 17387051

real	0m4.358s
user	0m4.688s
sys	0m1.795s


/Volumes/My_Book_2/raw_data_backup/Shanghai/55_RDM16_Col0_background_sRNA/55/result_primary_13:49:24_N=517$ time cat 55.fq | grep ATCGTATGCCGTCTTCTGCTTGA | wc -l
 4025722

/Volumes/My_Book_2/raw_data_backup/Shanghai/55_RDM16_Col0_background_sRNA/55/result_primary_13:52:14_N=518$ time cat 55.fq | grep TTCGTATGCCGTCTTCTGCTTGA | wc -l
 2567332

/Volumes/My_Book_2/raw_data_backup/Shanghai/55_RDM16_Col0_background_sRNA/55/result_primary_13:52:22_N=519$  cat 55.fq | grep CTCGTATGCCGTCTTCTGCTTGA | wc -l
 8125673

/Volumes/My_Book_2/raw_data_backup/Shanghai/55_RDM16_Col0_background_sRNA/55/result_primary_13:52:30_N=520$  cat 55.fq | grep GTCGTATGCCGTCTTCTGCTTGA | wc -l
 2588848


#| fastx_trimmer -f4 -l37 |  fastx_clipper -a CTGTAGGCAC -l 18 -c -v | fastq_to_fasta -v  -o  *.fa

#/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/small_RNA_data_from_SH_server/clean_data_12:22:23_N=520$  time zcat ../raw_fq/nrpd1.fq.gz | fastx_trimmer -f4 -l37 |  fastx_clipper -a CTGTAGGCAC -l 18 -c -v | fastq_to_fasta -v  -o nrpd1_clean_18_32.fa


adaptor: TCGTATGCCGTCTTCTGCT

time cat *.fq |  fastx_clipper -a TCGTATGCCGTCTTCTGCT -l 18 -c -v | fastq_to_fasta -v  -o  *.fa

time cat debug.fq |  fastx_clipper -a TCGTATGCCGTCTTCTGCT -l 18 -c -v | fastq_to_fasta -v  -o  debug_18-30.fa

##############
# 1,real trim
##############
time cat Col.fq  |  fastx_clipper -a TCGTATGCCGTCTTCTGCT -l 18 -c -v | fastq_to_fasta -v -o  Col0_WT_CF_clean_18_30.fa
time cat shh1.fq |  fastx_clipper -a TCGTATGCCGTCTTCTGCT -l 18 -c -v | fastq_to_fasta -v -o  shh1_dtf1_CF_clean_18_30.fa
time cat 55.fq   |  fastx_clipper -a TCGTATGCCGTCTTCTGCT -l 18 -c -v | fastq_to_fasta -v -o  55_rdm16-2_CF_clean_18_30.fa



# time bowtie -t -f -p 8 -v 0 -k 100   /Users/tang58/DataBase/TopHat/Ensembl/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BowtieIndex/genome *.fa  -S Col-0_trim_adpter.fq.trimmed.single.18-32.fa.bowtie_v0k100.sam



##############
# 2, length distribution
##############
# time cat *.fa | awk '{if(NR%2 == 0) cnt[length($1)] ++} END { for (x in cnt){print x,cnt[x]} }'  | sort -k1,1n | perl -lane 'print join("\t", @F)' > *txt #nrpd1_clean_18_32_length_stat.txt

 time cat Col0_WT_CF_clean_18_30.fa | awk '{if(NR%2 == 0) cnt[length($1)] ++} END { for (x in cnt){print x,cnt[x]} }'  | sort -k1,1n | perl -lane 'print join("\t", @F)' >  Col0_WT_CF_clean_18_30_length_stat.txt
 time cat shh1_dtf1_CF_clean_18_30.fa | awk '{if(NR%2 == 0) cnt[length($1)] ++} END { for (x in cnt){print x,cnt[x]} }'  | sort -k1,1n | perl -lane 'print join("\t", @F)' > shh1_dtf1_CF_clean_18_30_length_stat.txt
 time cat 55_rdm16-2_CF_clean_18_30.fa | awk '{if(NR%2 == 0) cnt[length($1)] ++} END { for (x in cnt){print x,cnt[x]} }'  | sort -k1,1n | perl -lane 'print join("\t", @F)' > 55_rdm16-2_CF_clean_18_30_length_stat.txt

/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background_16:56:38_N=578$ wc -l *fa
 34397476 55_rdm16-2_CF_clean_18_30.fa
 37610650 Col0_WT_CF_clean_18_30.fa
 27176946 shh1_dtf1_CF_clean_18_30.fa
 99185072 total


#########
# 3, bowtie in Terminal
###################


#########
# 4, split
###################
/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background_18:25:39_N=558$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.1.1.pl ln_unsorted_bam/dtf1_sRNA_CF_clean_bowtie_v0k100.bam splited_bam/ dtf1_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	13m40.199s
user	17m12.145s
sys	0m22.203s



/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_18:27:41_N=528$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.1.1.pl ../ln_unsorted_bam/rdm16-2_55_sRNA_CF_clean_bowtie_v0k100.bam . dm16-2_55_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	21m45.880s
user	27m18.926s
sys	0m35.000s


/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_18:27:45_N=540$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.1.1.pl ../ln_unsorted_bam/Col0_sRNA_CF_clean_bowtie_v0k100.bam . Col0_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	25m11.824s
user	31m27.023s
sys	0m40.312s



##################
# cal 500bp HNA db
###########################
time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl 
# new_HM_WT_v0.1_check_whole_no_structural_RNA_500bp_HNA_all.txt


/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_17:51:09_N=556$ 
time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl Col0_sRNA_CF_clean_bowtie_v0k100_without_structural_RNA.bam        WT_CF  15017949  WT_CF_v0.1_check_whole_no_structural_RNA_500bp_HNA_all.txt
time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl dtf1_sRNA_CF_clean_bowtie_v0k100_without_structural_RNA.bam        dtf1_CF  8687211  dtf1_CF_v0.1_check_whole_no_structural_RNA_500bp_HNA_all.txt
time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl rdm16-2_55_sRNA_CF_clean_bowtie_v0k100_without_structural_RNA.bam  rdm16-2_CF  13997762  rdm16-2_CF_v0.1_check_whole_no_structural_RNA_500bp_HNA_all.txt

######################################################################################################################################################################################################


#########
# 4, split to exclude overlapped sRNA
###################

/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_22:06:22_N=508$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.0.1.pl ../ln_unsorted_bam/Col0_sRNA_CF_clean_bowtie_v0k100.bam excluding_overlapping_stru_RNA/ Col0_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	24m8.003s
user	29m35.883s
sys	0m34.802s

/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_22:08:01_N=503$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.0.1.pl ../ln_unsorted_bam/rdm16-2_55_sRNA_CF_clean_bowtie_v0k100.bam excluding_overlapping_stru_RNA/ rdm16-2_55_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	20m43.064s
user	25m27.574s
sys	0m32.060s

/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam_22:07:36_N=504$ time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step1_split_bam2_twoBam_add_NH_tag_v0.0.1.pl ../ln_unsorted_bam/dtf1_sRNA_CF_clean_bowtie_v0k100.bam excluding_overlapping_stru_RNA/ dtf1_sRNA_CF_clean_bowtie_v0k100
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.
[samopen] SAM header is present: 7 sequences.

real	13m11.473s
user	15m55.795s
sys	0m20.162s

####################
# 5, cal 500bp HNA db
##################
time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl Col0_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam        WT_CF      14938655    WT_CF_noOver_structural_RNA_500bp_HNA_all.txt
real	15m42.418s
user	16m48.589s
sys	0m12.150s


time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl dtf1_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam        dtf1_CF    8634717     dtf1_CF_noOver_structural_RNA_500bp_HNA_all.txt
real	8m35.174s
user	9m11.876s
sys	0m6.709s

time perl ~/scripts_all/perl_code/Purdue/small_RNA_analysis/step2_cal_HNA_for_500bp_bin_v0.1.pl rdm16-2_55_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam  rdm16_2_CF 13920244    rdm16_2_CF_noOver_structural_RNA_500bp_HNA_all.txt
real	14m21.288s
user	15m22.885s
sys	0m11.232s


###############
# cal sRNA in DMR
#####################
~/misc/Chaofeng_Huang/20121007_rdm16/20130904_Col0_background_paper/sRNA_analysis/sRNA_in_DMR_16:00:20_N=504$

time perl /Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/step2.1_cal_HNA_for_non_overlaping_interval_v0.0_all_length.pl rdm16-2A_vs_Col0_CF_hypo_P0.05_reduced_boundary_Both_dep4_WinSize200_sliding50_gap100_iniCut2_repCut5_RDM16-2_detail_len100_CG30_CHG15_CHH10_annotated_TAIR10.txt  Col0_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam        WT_CF      14938655  rdm16_2A_vs_Col0_CF_hypo_915loci_sRNA_in_WTCF.txt


time perl /Users/tang58/scripts_all/perl_code/Purdue/small_RNA_analysis/step2.1_cal_HNA_for_non_overlaping_interval_v0.0_all_length.pl rdm16-2A_vs_Col0_CF_hypo_P0.05_reduced_boundary_Both_dep4_WinSize200_sliding50_gap100_iniCut2_repCut5_RDM16-2_detail_len100_CG30_CHG15_CHH10_annotated_TAIR10.txt  rdm16-2_55_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam  rdm16_2_CF 13920244   rdm16_2A_vs_Col0_CF_hypo_915loci_sRNA_in_rdm16_2CF.txt


###############
# sort bam file
#####################
/Volumes/My_Book/20130206_Huiming_smallRNA_Seq/Chaofeng_rdm16-2_Col0_background/splited_bam/excluding_overlapping_stru_RNA_10:08:50_N=502$ time samtools sort dtf1_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam dtf1_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA_sorted
[bam_sort_core] merging from 8 files...

real	3m15.902s
user	2m36.922s
sys	0m5.458s

 time samtools sort Col0_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam Col0_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA_sorted
[bam_sort_core] merging from 15 files...

real	6m10.279s
user	5m5.342s
sys	0m8.392s

time samtools sort rdm16-2_55_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA.bam rdm16_2_55_sRNA_CF_clean_bowtie_v0k100_NotOver_structural_RNA_sorted
[bam_sort_core] merging from 14 files...

real	6m5.424s
user	4m55.819s
sys	0m8.197s

# viswa's small RNA data
# (1) check miRNA expression level in Ha27 and Sta1 vs ColRD
# (2) check siRNA accumulation on luc gene
#first check quality
perl check_quality.pl
# from fastq to fasta (without removing any seqs)
perl change_fastq_to_fasta.pl
# remove adaptors
perl batch_remove_adaptor.pl
# calculate length distribution
perl length_dist.pl JKZ_ColRD_lane5/s_5_sequence_ge9.fa JKZ_Ha27_lane6/s_6_sequence_ge9.fa JKZ_Sta1_lane8/s_8_sequence_ge9.fa >clean_reads_length_dist.txt
# count number of miRNAs
 perl count_num_miRNAs.pl miRBase/release_17/ath_mature_T.fa JKZ_ColRD_lane5/s_5_sequence_ge9.fa JKZ_Ha27_lane6/s_6_sequence_ge9.fa  JKZ_Sta1_lane8/s_8_sequence_ge9.fa >num_miRNAs_in_libs.txt

# make a flat file so that each unique probe is included.
perl create_flat_file_ind_probe.pl ~/projects/tiling_array/bpmap/At35b_MR_v04-2_TAIR9_unique.tpmap >At35b_MR_TAIR10_ind_probe.flat

in R:
setwd("/media/disk/Heng_Zhang_tiling_array/2_13_2012/make_cdf")
source("flat2Cdf.R")
flat2Cdf(file="At35b_MR_TAIR10_ind_probe.flat", chipType="At35b_MR_v04",tags="ind_probe", splitn=6, rows=2560, cols=2560, verbose=8)

mv make_cdf/At35b_MR_v04,ind_probe.cdf annotationData/chipTypes/At35b_MR_v04/

in R:
setwd("/media/disk/Heng_Zhang_tiling_array/2_13_2012")
source("aroma_commands.R")
source("readCelFiles_normalized.R")
mv /media/disk/Heng_Zhang_tiling_array/2_13_2012/probeData/eight_samples,RBC,QN/At35b_MR_v04/*_RBC_QN.txt normalized_data

under normalized_data
perl ../extract_probe_vals.pl ~/projects/tiling_array/bpmap/At35b_MR_v04-2_TAIR9_unique.tpmap Col0_1.CEL_RBC_QN.txt 13028_1.CEL_RBC_QN.txt 24055_1.CEL_RBC_QN.txt chr40_1.CEL_RBC_QN.txt nrpd1a_1.CEL_RBC_QN.txt nrpd1b_1.CEL_RBC_QN.txt pkl_1_1.CEL_RBC_QN.txt shh1_1_1.CEL_RBC_QN.txt >../normalized_probe_intensities.txt

# find differentially expressed regions
perl compare_two_samples.pl normalized_probe_intensities.txt 2 3 0 >Col0_vs_13028_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 3 1 >Col0_vs_13028_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 4 1 >Col0_vs_24055_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 4 0 >Col0_vs_24055_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 5 0 >Col0_vs_chr40_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 5 1 >Col0_vs_chr40_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 6 1 >Col0_vs_nrpd1a_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 6 0 >Col0_vs_nrpd1a_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 7 0 >Col0_vs_nrpd1b_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 7 1 >Col0_vs_nrpd1b_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 8 1 >Col0_vs_pkl_down.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 8 0 >Col0_vs_pkl_up.txt
perl compare_two_samples.pl normalized_probe_intensities.txt 2 9 0 >Col0_vs_shh1_up.txt
rar a diff_expr_regions.rar *_down.txt
rar a diff_expr_regions.rar *_up.txt

#R --slave --vanilla --args indir list_pre binSize half_bin_num  outpre
# < /Users/tang58/scripts_all/perl_code/Purdue/Histone_src/histone_score/Nov21_read_density_RPKM/draw_DMRs_histone_RPKM_one_fig.r

args = commandArgs( trailingOnly =  T)  
indir		= args[1]
list_pre	= args[2]
binSize		= args[3]; binSize = as.numeric(binSize)
half_bin_num= args[4]; half_bin_num = as.numeric(half_bin_num)
png_pre		= args[5]


###########
# set up tag number for RPKM calculation
#############
library(hash)
tag_num_file = "/Volumes/Macintosh_HD_2/Kai_project/cross_talk_of_RdDM_and_methylation/Jacobsen_PNAS_method_v2/downstream/connect_with_histone/IDM1_paper_list/ReadsN_as_PLOS_method/bam_tag_number.txt"
tag_num = read.table(tag_num_file, head=T)
h = hash(tag_num[,1], tag_num[,2])
# h[[  ]]

##########
# read files and get histome lable
################
setwd(indir)
files = dir(pattern= list_pre )

file_num = length(files)
list_data = list()
for (i in 1:file_num){ 
	list_data[[i]] = read.delim(files[i], head=T, sep = "\t") 
}
bam_lables = sub(".+_in_(.+)_half.+.txt", "\\1", files, perl =T)
#ros1_4_hyper_P0.01_reduced_boundary_both_depth100_WinSize1000_gap1000_initialCutoff5_reportCutoff10_in_H3K18Ac_halfBinN25_bin200bp_ReadsN.txt
#print (bam_lables)

##########
# color
#########
#(my_col = palette(rainbow(11)))
#my_col = rainbow(11)
# "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow"  "gray" "purple"
#"#FF8B00" TuHuang
#"#808000" Olive
##008080 teal
my_col = c(  "red",     "green3",  "blue",    "cyan",    "magenta", "#008080",  "gray", "purple", "#FF8B00" , "#808000","black")
cex_lab = 2
cex_axis = 1.8
cex_axis_plot = 1.5
mgp =  c(5, 0.8, 0)

width_val = 9# 3.8; 
height_val = 7#4;
units_val = "in"; res_val = 500; pointsize_val =8;

#setwd(outdir)
#file_name = paste(outdir, "/", png_pre , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")
file_name = paste(png_pre , "_w", width_val, "_h", height_val, "_res", res_val, "_pt", pointsize_val, ".png", sep="")
png(file_name, width = width_val, height = height_val, units = units_val , res = res_val, pointsize = pointsize_val );

par(oma = c(10, 3, 1, 0) ) #down left up rigth
par(mar = c(2, 10 , 2, 1) + 0.1)


lwd_abline = 1.5
lwd_val   = 2#1.2
las_val = 1; # axis labels horizontal or vertical


########
# cal RPKM and record max
########
#RPKM = (10^9 * C)/(N * L)

max_all = 0
#list_RPKM = list()
list_RPKM_frame = data.frame()
list_RPKM_mean = list()
for (i in 1:file_num){ 
	label =  bam_lables[i]
	mapped_reads = h[[label]]
	list_RPKM_frame = (list_data[[i]][,3:(2*half_bin_num + 3)]) * 10^9 / (binSize * mapped_reads)
	list_RPKM_mean[[i]] = colMeans(list_RPKM_frame)
	max_tmp = max(list_RPKM_mean[[i]] )
	if( max_tmp > max_all ){
		max_all = max_tmp
	}	
}

######
mid_point = half_bin_num+ 1# 101#41

i = 1

plot(	
		#x[,1],x[,i]/x[,12], 
		1:( 2*half_bin_num + 1 ),list_RPKM_mean[[i]],
		ylim = c(0, max_all),
		type="l",#b,o,l
		col= my_col[i], #"blue", 
		xaxt="n", las = las_val,xlab="",
		ylab = "ChIP-seq read density (RPKM)", #colnames(x)[i], 
		cex.lab= cex_lab, cex.axis = cex_axis_plot, mgp = mgp,
		lwd = lwd_val
)

abline(v=mid_point, lty= "dotted") 

axis(1, at = c(1,mid_point,2*half_bin_num + 1), labels=c("-5kb","Midpoint","+5kb"), cex.axis=cex_axis)

for (i in 2:file_num){
	lines(1:( 2*half_bin_num + 1 ),list_RPKM_mean[[i]], col = my_col[i], lwd = lwd_val )
}

par(xpd=NA)
legend (   13, -2,
legend= bam_lables, 
col = my_col,
lwd = 2,
bty = "n",
cex = 1.5,
ncol = 3
)

dev.off()


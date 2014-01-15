
setwd("/media/disk/Heng_Zhang_tiling_array/2_13_2012/probeData/eight_samples,RBC,QN/At35b_MR_v04")
library(affy)
library(affxparser)
### read .cel files
celFiles <- list.files(pattern="[.](c|C)(e|E)(l|L)$")
for (i in 1:length(celFiles)) {
	celData <- readCel(celFiles[i], readHeader=TRUE, readXY=TRUE)
    data_table <- cbind(celData$x, celData$y, celData$intensities)
	outfile <- paste(celFiles[i], "_RBC_QN.txt", sep="")
	write.table(data_table, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}
system("mv /media/disk/Heng_Zhang_tiling_array/2_13_2012/probeData/eight_samples,RBC,QN/At35b_MR_v04/*_RBC_QN.txt /media/disk/Heng_Zhang_tiling_array/2_13_2012")
# then change file names


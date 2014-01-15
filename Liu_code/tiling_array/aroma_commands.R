library(aroma.affymetrix)
setwd("/media/disk/Heng_Zhang_tiling_array/2_13_2012")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
geneCdf <- AffymetrixCdfFile$byChipType("At35b_MR_v04,ind_probe")
chipType <- "At35b_MR_v04"
cs <- AffymetrixCelSet$byName("eight_samples",cdf=geneCdf, chipType=chipType)
#check raw images.files located at reports/<data set name>/<data set tags>/<chipType>/spatial/<fullname>,<tags>.png

#lapply(cs, writeImage, tags="rainbow", palette=rainbow(256))

#quality assessment
chipNames <- cs$Names
pdf("rawdata_density_plot_pm_eight_samples.pdf")
plotDensity(cs, ylim=c(0,1.5), col=1:5, lwd=1, types="pm")
legend("topright", chipNames, col=1:5, lwd=1)
dev.off()

# backgroud correction
bc <- RmaBackgroundCorrection(cs)

csBC <- process(bc, verbose=verbose)

# quantile normalization
qn <- QuantileNormalization(csBC, typesToUpdate="pm")

csN <- process(qn, verbose=verbose)



# ###maser for rmats
# ### https://www.bioconductor.org/packages/devel/bioc/vignettes/maser/inst/doc/Introduction.html#1_Overview_of_maser
# Maser Tutorial online
#https://bioconductor.org/packages/release/bioc/vignettes/maser/inst/doc/Introduction.html#2_Importing_rMATS_events
BiocManager::install("maser")

library(maser)
library(rtracklayer)

#Path to 468_0hours vs 468_12hours
path <- "./468_0_v_468_12"

#Calling the Rmats file we are using
Met <- maser(path, c("MB468_0mins", "MB468_720mins"), ftype = "JCEC")
Met

#head of summary of the rmats output
head(summary(Met, type = "SE")[, 1:8])

Met_filt <- filterByCoverage(Met, avg_reads = 5)
#Low coverage splicing junctions are commonly found in RNA-seq data and lead to low confidence PSI levels. We can remove low coverage events using filterByCoverage(), which may signficantly reduced the number of splicing events.


Met_top <- topEvents(Met_filt, fdr = 0.05, deltaPSI = 0.1)
Met_top
#The function topEvents() allows to select statistically significant events given a FDR cutoff and minimum PSI change. Default values are fdr = 0.05 and deltaPSI = 0.1 (ie. 10% minimum change).

#An overview of significant events can be obtained using either dotplot() or volcano() functions, specifying FDR levels, minimum change in PSI between conditions and splicing type. Significant events in each condition will be highlighted.
volcano(Met, fdr = 0.05, deltaPSI = 0.1, type = "SE")
examplea <- volcano(Met_top, fdr = 0.05, deltaPSI = 0.1, type = "SE")
examplea

#theme() allows you to change aspects of the volcano plot
examplea + theme(axis.title.x = element_text(size=18)) + theme(axis.title.y = element_text(size=18))+ theme(legend.key=element_rect(fill='bisque')) + guides(colour = guide_legend(override.aes = list(size=4))) + theme(legend.text=element_text(size=15)) + theme(axis.title = element_text(size = 18)) + 
theme(legend.title=element_blank()) + theme(axis.text.x  = element_text(size=16)) + ggtitle("MB468_0mins vs MB468_720mins (Skipped Exon)") + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(axis.text.y  = element_text(size=16)) + theme(plot.title = element_text(size=18)) +
xlim(-0.75, 0.75)

#If only significant events should be plotted, then use topEvents() combined with volcano() or dotplot() for visualization.
example1 <- volcano(Met_top, fdr = 0.05, deltaPSI = 0.1, type = "RI")
example1

example1 + theme(axis.title.x = element_text(size=18)) + theme(axis.title.y = element_text(size=18))+ theme(legend.key=element_rect(fill='bisque')) + guides(colour = guide_legend(override.aes = list(size=4))) + theme(legend.text=element_text(size=15)) + theme(axis.title = element_text(size = 18)) + 
  theme(legend.title=element_blank()) + theme(axis.text.x  = element_text(size=16)) + ggtitle("MB468_0 vs MB468_720 Intron Retention") + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(axis.text.y  = element_text(size=16)) + theme(plot.title = element_text(size=18))

# type = refers to different splicing events
volcano(Met_filt, fdr = 0.05, deltaPSI = 0.1, type = "RI")
volcano(Met_filt, fdr = 0.05, deltaPSI = 0.1, type = "MXE")
volcano(Met_filt, fdr = 0.05, deltaPSI = 0.1, type = "A5SS")
volcano(Met_filt, fdr = 0.05, deltaPSI = 0.1, type = "A3SS")

#The breakdown of splicing types can be plotted using splicingDistribution() and desired significance thresholds. Please refer to help pages for examples on how to use these functions.
splicingDistribution(Met_filt, fdr = 0.05, deltaPSI = 0.1)

#dotplot visualization
dotplot(Met_top, type = "SE")
####


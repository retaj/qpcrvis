library(qpcrvis)

# load demo data
xls.files <- dir(system.file("extdata/", package="qpcrvis"), full.names = T)

# this will load all xls.files into a list
pcrs <- lapply(xls.files, readPCR)

# OK, we have 2 plates with 4 samples and 8 genes each (reporting correct number of wells is currently too much for me)
pcrs

# we can have a look at them separately
plotRQ(pcrs[[1]])
plotRQ(pcrs[[2]])

# naming is also cool if you prefer to reference stuff using $ notation
names(pcrs) <- c("my_first_qPCR", "my_second_qPCR")

# now, as both qPCR runs are a part of one experiment, let's try to merge the plates and redo the stats
pcrs.merged <- mergePCR(list(pcrs$my_first_qPCR@data, pcrs$my_second_qPCR@data), "WT1", "geneA")
# not quite what we expected, people are horrible and use inconsistent gene and sample naming styles
pcrs.merged
# also visible from the summarized list
qPCRsummary(pcrs)

# let's rename the samples on the second plate...
pcrs$my_second_qPCR <- renameSamples(pcr = pcrs$my_second_qPCR, old = levels(pcrs$my_second_qPCR@data$sample), new = c("KO3", "KO4", "WT3", "WT4"))
# ...genes too
pcrs$my_second_qPCR <- renameTargets(pcr = pcrs$my_second_qPCR, old = levels(pcrs$my_second_qPCR@data$target), new = paste0("gene", LETTERS[1:8]))
# looks better
qPCRsummary(pcrs)
pcrs.merged <- mergePCR(list(pcrs$my_first_qPCR@data, pcrs$my_second_qPCR@data), "WT1", "geneA")
pcrs.merged

# kk, plot
plotRQ(pcrs.merged)

# sample ordering is shitty, improve
pcrs.merged.ord <- reorderSamples(pcrs.merged,
                                  old = levels(pcrs.merged@data$sample),
                                  new = levels(pcrs.merged@data$sample)[c(grep("WT", levels(pcrs.merged@data$sample)), grep("KO", levels(pcrs.merged@data$sample)))]
)

plotRQ(pcrs.merged.ord)
plotRQ(pcrs.merged.ord, plottype = "point")

# R has no idea what wild type and knockout are...but we can set the design manually
pcrs.merged.ord <- setDesign(pcrs.merged.ord, groups = rep(c("WT", "KO"), each = 4))
plotRQ(pcrs.merged.ord)
plotRQ(pcrs.merged.ord, plottype = "point")

# want to reorder genes? wait (or hack) until it's implemented
# dropping genes is fine tho
pcrs.merged.ord <- keepTargets(pcr = pcrs.merged.ord, keep = c("geneC", "geneD", "geneE", "geneF", "geneG"))
plotRQ(pcrs.merged.ord)

# just as dropping samples, but will break the color scheme, so design shall be redefined (will fix)
pcrs.merged.ord <- keepSamples(pcr = pcrs.merged.ord, keep = c("WT1", "WT3", "KO3", "KO4"))
pcrs.merged.ord <- setDesign(pcrs.merged.ord, groups = rep(c("WT", "KO"), each = 2))
plotRQ(pcrs.merged.ord)

# don't like the default plot appearance? get raw ggplot2 objects and edit
p <- plotRQ(pcrs.merged.ord, minimal = TRUE)
print(p)

library(ggplot2)
p <- p + ggtitle("i like defaults and funky y axes") + ylim(0,100)
print(p)

# this plots mean Cts across samples, but needs to be improved.
plotRawCt(pcrs$my_first_qPCR, raw = FALSE)

# ---------------------------------------------------------------------------- #
#' setDesign
#'
#' set design of a qPCR experiment
#'
#' details
#'
#' @param pcr qPCR object to work on
#' @param groups character vector of sample groups
#'
#' @importFrom xlsx read.xlsx
#' @importFrom data.table data.table setnames set
#' @export
setGeneric(
  name="setDesign",
  def=function(pcr, ...) {
    standardGeneric("setDesign")
  }
)
setMethod("setDesign",
          signature("qPCR"),
          definition=function(pcr, groups) {
            if (length(groups) != length(unique(pcr@data$sample))) {
              stop("cannot assign ", length(groups), " group label(s) to ", length(unique(pcr@data$sample)), " samples.")
            }
            pcr@design <- factor(groups, levels=unique(groups))

            return(pcr)
          }
)


# ---------------------------------------------------------------------------- #
#' renameSamples
#'
#' set new sample names for a qPCR object
#'
#' @param pcr qPCR object to work on
#' @param groups character vector of sample groups
#'
#' @importFrom xlsx read.xlsx
#' @importFrom data.table data.table setnames set
#' @export
setGeneric(
  name="renameSamples",
  def=function(pcr, old, new) {
    standardGeneric("renameSamples")
  }
)
setMethod("renameSamples",
          signature("qPCR"),
          definition=function(pcr, old, new) {

            # stay on the character side of life
            old <- as.character(old)

            # exceptions!
            if (length(old) != length(new))
              stop("different number of old and new sample names supplied: ", length(old), " vs. ", length(new), ".")

            if (!identical(sort(old), sort(as.character(levels(pcr@data$sample)))))
              stop("old sample names not identical to sample names in qPCR object.")

            # reorder first
            pcr@data$sample <- factor(pcr@data$sample, levels=old)
            # aaand rename
            #pcr@data$sample <- factor(pcr@data$sample, levels=new)
            levels(pcr@data$sample) <- new

            return(pcr)
          }
)





# ---------------------------------------------------------------------------- #
#' relExp
#'
#' calculate relative expression for a qPCR run
#'
#' details
#'
#' @param pcr qPCR object to work on
#' @param groups character vector of sample groups
#'
#' @importFrom xlsx read.xlsx
#' @importFrom data.table data.table setnames
#' @importFrom reshape2 melt dcast
#' @export
setGeneric(
  name="relExp",
  def=function(pcr, ...) {
    standardGeneric("relExp")
  }
)
setMethod("relExp",
          signature("qPCR"),
          definition=function(pcr, ref_sample, ref_target) {

            if (!(ref_sample %in% pcr@raw.data$sample) | !(ref_target %in% pcr@raw.data$target)) {
              stop("reference sample or target gene not present in raw data")
            }

            DT.raw <- pcr@raw.data[,.(sample, target, Ct)]
            DT <- cbind(data.table(melt(dcast(DT.raw, target~sample, mean))),
                        data.table(melt(dcast(DT.raw, target~sample, sd))))
            setnames(DT, c("target", "sample", "Ct_mean", "target2", "sample2", "Ct_sd"))
            if (identical(DT$target, DT$target2) & identical(DT$sample, DT$sample2)) {
              DT <- DT[,.(sample, target, Ct_mean, Ct_sd)]
            } else {
              stop("something is fishy with sample and/or target labels. blame the developer.")
            }
            DT <- DT[order(target)]

            # the table has to be properly sorted for this!
            DT$dCt <- DT$Ct_mean -  DT[target == ref_target]$Ct_mean
            DT$dCt_sd <- sqrt(DT$Ct_sd**2 + DT[target == ref_target]$Ct_sd**2)

            ddcts <- numeric()
            for (i in unique(DT$target)) {
              ddcts <- append(ddcts, DT[target==i]$dCt - DT[target==i & sample==ref_sample]$dCt)
            }
            DT$ddCt <- ddcts

            DT$RQ <- 2**(-DT$ddCt)
            DT$RQmin <- 2**(-DT$ddCt - DT$dCt_sd)
            DT$RQmax <- 2**(-DT$ddCt + DT$dCt_sd)

            pcr@data <- DT

            return(pcr)
          }
)








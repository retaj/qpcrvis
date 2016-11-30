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

            if (!(ref_sample %in% pcr@data$sample) | !(ref_target %in% pcr@data$target)) {
              stop("reference sample or target gene not present in raw data")
            }

            DT.raw <- pcr@data[,.(well, sample, target, Ct)]
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

            pcr@pdata <- DT

            return(pcr)
          }
)








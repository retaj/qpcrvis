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
#' renameTargets
#'
#' set new target names for a qPCR object
#'
#' @param pcr qPCR object to work on
#' @param old character vector of target names to be replaced
#' @param new character vector of new target names
#'
#' @importFrom data.table data.table setnames set
#' @export
setGeneric(
  name="renameTargets",
  def=function(pcr, old, new) {
    standardGeneric("renameTargets")
  }
)
setMethod("renameTargets",
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
            pcr@data$target <- factor(pcr@data$target, levels=old)
            # aaand rename
            #pcr@data$sample <- factor(pcr@data$sample, levels=new)
            levels(pcr@data$target) <- new

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

            # Ct values can be NA, I can still drop them and calculate mean and SD from the rest
            if (nrow(DT.raw[is.na(Ct)]) > 0) {
              warning("raw Ct NA values will be ignored:\n", paste(capture.output(print(DT.raw[is.na(Ct)])), collapse="\n"))
              DT.raw <- DT.raw[!is.na(Ct)]
            }

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



# ---------------------------------------------------------------------------- #
#' mergePCR
#'
#' merge 2 qPCR objects
#'
#' details
#'
#' @param pcr1,pcr2  qPCR objects to work on
#' @param ref_sample sample to normalize to
#'
#' @export
setGeneric(
  name="mergePCR",
  def=function(pcrs, ref_sample, ref_gene) {
    standardGeneric("mergePCR")
  }
)
setMethod("mergePCR",
          signature("list"),
          definition=function(pcrs, ref_sample, ref_gene) {
            # TODO: write this using ... later

            if (!hasArg(ref_sample))
              stop("please provide a reference sample!")
            if (!hasArg(ref_gene))
              stop("please provide a gene used to normalize plates!")

            #pcrs <- lapply(list(...), function(x) x@data)
            # TODO: check if both are normalized to the same ref_target!

            #DT.merged <- rbind(pcr1@data, pcr2@data)
            DT.merged <- do.call("rbind", pcrs)
            levels(DT.merged$target) <- tolower(levels(DT.merged$target))

            # dCts and dCt_sds are fine, i need to
            # recalculate ddCt and RQ/RQmin/RQmax
            DT <- DT.merged[,.(sample, target, Ct_mean, Ct_sd, dCt, dCt_sd)][order(target)]

            ddcts <- numeric()
            for (i in unique(DT$target)) {
              ddcts <- append(ddcts, DT[target==i]$dCt - DT[target==i & sample==ref_sample]$dCt)
            }
            DT$ddCt <- ddcts

            DT$RQ <- 2**(-DT$ddCt)
            DT$RQmin <- 2**(-DT$ddCt - DT$dCt_sd)
            DT$RQmax <- 2**(-DT$ddCt + DT$dCt_sd)

            pcr <- new("qPCR", raw.data = data.table(),
                               metadata = data.table(X1='Endogenous Control',
                                                     #X2=names(which(table(DT[RQ==1]$target) == length(levels(DT$sample))))), # TODO: this is dirty af
                                                     X2=ref_gene), # TODO: this is dirty af
                               data     = DT
            )

            return(pcr)
          }
)

# ---------------------------------------------------------------------------- #
#' reorderSamples
#'
#' set the order of samples
#'
#' details
#'
#' @param pcr qPCR objects to work on
#' @param old old sample order
#' @param new new sample order
#'
#' @export
setGeneric(
  name="reorderSamples",
  def=function(pcr, old, new) {
    standardGeneric("reorderSamples")
  }
)
setMethod("reorderSamples",
          signature("qPCR"),
          definition=function(pcr, old, new) {
            # stay on the character side of life
            old <- as.character(old)

            # exceptions! TODO: extract this to a separate function to be used by all renames, reorders etc
            if (length(old) != length(new))
              stop("different number of old and new sample names supplied: ", length(old), " vs. ", length(new), ".")

            if (!identical(sort(old), sort(as.character(levels(pcr@data$sample)))))
              stop("old sample names not identical to sample names in qPCR object.")

            if (!(sum(old %in% new) == length(old)))
              stop("different samples in old and nww")

            # set new factors
            pcr@data$sample <- factor(as.character(pcr@data$sample), levels=new)

            return(pcr)
          }
)

# ---------------------------------------------------------------------------- #
#' keepTargets
#'
#' keep only selected targets
#'
#' details
#'
#' @param pcr  qPCR objects to work on
#' @param keep character vector of targets to keep
#'
#' @export
setGeneric(
  name="keepTargets",
  def=function(pcr, keep) {
    standardGeneric("keepTargets")
  }
)
setMethod("keepTargets",
          signature("qPCR"),
          definition=function(pcr, keep) {
            # stay on the character side of life
            keep <- unique(as.character(keep))

            # exceptions! TODO: extract this to a separate function to be used by all renames, reorders etc
            if (sum(!(keep %in% levels(pcr@data$target))) > 0)
              stop("unknown levels: ", paste0(keep[!(keep %in% levels(pcr@data$target))], collapse=", "), ".")

            # set new factors
            pcr@data <- pcr@data[target %in% keep]
            pcr@data$target <- factor(pcr@data$target, levels=levels(pcr@data$target)[levels(pcr@data$target) %in% keep])

            return(pcr)
          }
)


# ---------------------------------------------------------------------------- #
#' keepSamples
#'
#' keep only selected samples
#'
#' details
#'
#' @param pcr  qPCR objects to work on
#' @param keep character vector of samples to keep
#'
#' @export
setGeneric(
  name="keepSamples",
  def=function(pcr, keep) {
    standardGeneric("keepSamples")
  }
)
setMethod("keepSamples",
          signature("qPCR"),
          definition=function(pcr, keep) {
            # stay on the character side of life
            keep <- unique(as.character(keep))

            # exceptions! TODO: extract this to a separate function to be used by all renames, reorders etc
            if (sum(!(keep %in% levels(pcr@data$sample))) > 0)
              stop("unknown levels: ", paste0(keep[!(keep %in% levels(pcr@data$sample))], collapse=", "), ".")

            # set new factors
            pcr@data <- pcr@data[sample %in% keep]
            pcr@data$sample <- factor(pcr@data$sample, levels=levels(pcr@data$sample)[levels(pcr@data$sample) %in% keep])

            return(pcr)
          }
)

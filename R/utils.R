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

# ---------------------------------------------------------------------------- #
#' readPCR
#'
#' plot relative quantification as returned by the machine
#'
#' details
#'
#' @param filename
#'
#' @importFrom xlsx read.xlsx
#' @importFrom data.table data.table setnames
#' @export
setGeneric(
  name="plotRQ",
  def=function(pcr) {
    standardGeneric("plotRQ")
  }
)
setMethod("plotRQ",
          signature("qPCR"),
          definition=function(pcr) {
            print("plop")
          }
)

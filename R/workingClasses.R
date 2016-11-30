#' An S4 class to represent a qPCR experiment.
#'
#' @slot data     data.table with qPCR data
#' @slot metadata data.table with experimental metadata
setClass("qPCR",
         slots = c(data="data.table", metadata="data.table", design="factor")
)

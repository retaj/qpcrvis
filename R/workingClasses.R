# dirty patch to be able to use data.table in class slots
setOldClass(c('data.frame'))
setOldClass(c('data.table', 'data.frame'))

#' An S4 class to represent a qPCR experiment.
#'
#' @slot data     data.table with qPCR data
#' @slot metadata data.table with experimental metadata
#' @slot design   factor vector of sample groupings for barplot coloring
#' @slot pdata    processed data (e.g. recalculated relative expression)
#' @importFrom data.table data.table
#'
#' @name   qPCR-class
#' @rdname qPCR-class
#'
#' @exportClass qPCR
setClass("qPCR",
         slots = c( data      = "data.table",
                    metadata  = "data.table",
                    design    = "factor",
                    pdata     = "data.table")
)

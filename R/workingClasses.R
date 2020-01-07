# dirty patch to be able to use data.table in class slots
setOldClass(c('data.frame'))
setOldClass(c('data.table', 'data.frame'))

#' An S4 class to represent a qPCR experiment.
#'
#' @slot raw.data data.table with raw qPCR data
#' @slot metadata data.table with experimental metadata
#' @slot design   factor vector of sample groupings for barplot coloring
#' @slot data    processed data (e.g. recalculated relative expression)
#' @importFrom data.table data.table
#'
#' @name   qPCR-class
#' @rdname qPCR-class
#'
#' @exportClass qPCR
setClass("qPCR",
         slots = c( raw.data   = "data.table",
                    metadata   = "data.table",
                    design     = "factor",
                    ref_sample = "character",
                    ref_target = "character",
                    data       = "data.table")
)



# show qPCR
setMethod("show", "qPCR",
          function(object){
            #n.wells <- nrow(object@raw.data)
            #samples <- levels(object@raw.data$sample)
            #targets <- levels(object@raw.data$target)
            n.wells <- nrow(object@data)
            samples <- levels(object@data$sample)
            targets <- levels(object@data$target)

            summ = paste0(n.wells, "-well qPCR for ", length(targets), " targets in ", length(samples), " samples.", "\n",
                          "samples: ", paste(samples, collapse=", "), "\n",
                          "targets: ", paste(targets, collapse=", "), "\n",
                          "reference sample: ", object@ref_sample, "\n",
                          "reference target: ", object@ref_target
                          )
            message(summ)
          })

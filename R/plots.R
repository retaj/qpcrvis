#' @include workingClasses.R
NULL

#' ---------------------------------------------------------------------------- #
#' plotRQ
#'
#' plot relative quantification as returned by the machine
#'
#' details
#'
#' @param filename
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
setGeneric("plotRQ",
            function(pcr, raw=FALSE, minimal=FALSE, ...) standardGeneric("plotRQ"))
#' @aliases plotRQ,plotRQ-method
#' @rdname plotRQ
setMethod("plotRQ",
          signature("qPCR"),
          function(pcr, ...) {

            # munge data
            if (raw == FALSE) {
              DT <- pcr@data[target!=as.character(pcr@metadata$X2[pcr@metadata$X1=="Endogenous Control"])]
              DT <- unique(DT[,.(sample, target, RQ, RQmin, RQmax)])
            } else {
              DT <- pcr@raw.data[target!=as.character(pcr@metadata$X2[pcr@metadata$X1=="Endogenous Control"])]
              DT <- unique(DT[,.(sample, target, RQ, RQmin, RQmax)])
            }
            #DT <- DT[order(sample)]

            errs <- aes(ymax = RQmax, ymin=RQmin)
            dodge=position_dodge(width=0.75)

            p <- ggplot(DT, aes(x=target, y=RQ, fill=sample)) +
                  geom_bar(stat="identity", position=dodge, width=0.65) +
                  geom_errorbar(errs, position=dodge, width=0.33)

            if (minimal==FALSE) {
              p <- p +  scale_y_continuous(expand=c(0, 0)) +
                        xlab("target genes") +
                        ylab("relative expression") +
                        theme_classic() +
                        theme(axis.line.y  = element_line(size=.5),
                              axis.line.x  = element_line(size=.5),
                              axis.text.x  = element_text(size=18, family="Helvetica"),
                              axis.text.y  = element_text(size=18, family="Helvetica"),
                              axis.title.x = element_text(size=20, family="Helvetica"),
                              axis.title.y = element_text(size=20, family="Helvetica"),
                              plot.title   = element_text(size=24, face="bold", family="Helvetica"),
                              aspect.ratio = 0.5)

              if (length(pcr@design) > 0) {
                if (length(levels(pcr@design)) <= 3) {

                  # design is defined, color appropriately
                  palettes <- c("Blues", "Greens", "Purples")
                  design.table <- table(pcr@design)
                  color.values <- character()
                  for (i in 1:length(design.table)) color.values <- c(color.values, rev(brewer.pal(7, palettes[i]))[seq(2, 2+design.table[i], by=2)])
                  p <- p + scale_fill_manual(values=color.values)

                } else {

                  warning("only up to three pallettes are currently supported. set space_fill_manual manually or open github issues until I implement it :)")

                }
              }
            }


            return(p)
          }
)

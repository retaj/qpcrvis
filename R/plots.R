# ---------------------------------------------------------------------------- #
#' readPCR
#'
#' plot relative quantification as returned by the machine
#'
#' details
#'
#' @param filename
#'
#' @import ggplot2
#' @export
setGeneric(
  name="plotRQ",
  def=function(pcr, ...) {
    standardGeneric("plotRQ")
  }
)
setMethod("plotRQ",
          signature("qPCR"),
          definition=function(pcr, minimal=FALSE) {

            # munge data
            DT <- pcr@data[target!=as.character(pcr@metadata$X2[pcr@metadata$X1=="Endogenous Control"])]
            DT <- unique(DT[,.(sample, target, RQ, RQmin, RQmax)])
            #DT <- DT[order(sample)]

            errs <- aes(ymax = RQmax, ymin=RQmin)
            dodge=position_dodge(width=0.75)

            p <- ggplot(DT, aes(x=target, y=RQ, fill=sample)) +
                  geom_bar(stat="identity", position=dodge, width=0.65) +
                  geom_errorbar(errs, position=dodge, width=0.33)

            if (minimal==FALSE) {
              p <- p +  xlab("target genes") +
                        ylab("relative expression") +
                        theme_classic() +
                        theme(axis.line.y=element_line(size=.5),
                          axis.line.x=element_line(size=.5),
                          axis.text.x=element_text(size=18,  family="Helvetica"),
                          axis.text.y=element_text(size=18,  family="Helvetica"),
                          axis.title.x=element_text(size=20, family="Helvetica"),
                          axis.title.y=element_text(size=20, family="Helvetica"),
                          aspect.ratio=0.5,
                          plot.title=element_text(size=24, face="bold", family="Helvetica"))
            }

            return(p)
          }
)

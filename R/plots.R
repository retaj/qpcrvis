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
#' @importFrom data.table data.table
#' @export
setGeneric("plotRQ",
            function(pcr, raw=FALSE, minimal=FALSE, plottype=c("bar", "point"), reverseGroups=FALSE, ...) standardGeneric("plotRQ"))
#' @aliases plotRQ,plotRQ-method
#' @rdname plotRQ
setMethod("plotRQ",
          signature("qPCR"),
          function(pcr, ...) {

            plottype <- match.arg(plottype)

            if (reverseGroups == TRUE) {
              warning("reverseGroups suspended until further notice.")
            }
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

            if (plottype=="point") {
              # TODO: write this robustly! (fetch values from design)
              DT$group <- factor(ifelse(grepl("ko", DT$sample, ignore.case=TRUE), "KO", "WT"), levels=c("WT", "KO"))

              meanz <- cbind(rbind(data.table(aggregate(RQ ~ target, data=DT[group=="WT"], FUN=mean)),
                                   data.table(aggregate(RQ ~ target, data=DT[group=="KO"], FUN=mean))),
                             group=factor(rep(c("WT", "KO"), each=nrow(DT)/length(levels(DT$sample))), levels=c("WT", "KO")))

              # if (reverseGroups == TRUE) {
              #   DT$group    <- factor(DT$group,    levels=rev(levels(DT$group)))
              #   meanz$group <- factor(meanz$group, levels=rev(levels(meanz$group)))
              # }

              dodge <- position_jitterdodge(jitter.width = .2, dodge.width = .75)

              p <- ggplot(DT, aes(x=target, y=RQ, group=group, color=group)) +
                    geom_pointrange(errs, position=dodge) +
                    geom_crossbar(aes(ymin = RQ, ymax = RQ), data=meanz, position=position_dodge(width=0.75), width=0.5)

            } else {
              dodge=position_dodge(width=0.75)

              p <- ggplot(DT, aes(x=target, y=RQ, fill=sample)) +
                geom_bar(stat="identity", position=dodge, width=0.65) +
                geom_errorbar(errs, position=dodge, width=0.33)
            }

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
                              plot.title   = element_text(size=24, family="Helvetica", face="bold", hjust=0.5),
                              legend.text  = element_text(size=16, family="Helvetica"),
                              legend.title = element_blank(),
                              aspect.ratio = 0.5)

              if (length(pcr@design) > 0) {
                if (length(levels(pcr@design)) <= 3) {

                  # design is defined, color appropriately
                  palettes <- c("Greens", "Blues", "Purples")
                  design.table <- table(pcr@design)
                  color.values <- character()
                  for (i in 1:length(design.table)) color.values <- c(color.values, rev(brewer.pal(9, palettes[i]))[1:design.table[i]])
                  if (plottype=="point") {
                    # TODO: dirty override
                    p <- p + scale_color_manual(values=rev(c("#2171B5", "#238B45")))
                  } else {
                    p <- p + scale_fill_manual(values=color.values)
                  }

                } else {

                  warning("only up to three pallettes are currently supported. set space_fill_manual manually or open github issues until I implement it :)")

                }
              }
            }


            return(p)
          }
)


#' ---------------------------------------------------------------------------- #
#' plotRawCt
#'
#' plot raw Ct values for reference
#'
#' details
#'
#' @param pcr qPCR object to work on
#'
#' @import ggplot2
#' @export
setGeneric("plotRawCt",
           function(pcr, raw=TRUE) standardGeneric("plotRawCt"))
setMethod("plotRawCt",
          signature("qPCR"),
          function(pcr, raw=TRUE) {

            if (raw == TRUE) {
              DT <- pcr@raw.data[,.(sample, target, Ct)]
              p <- ggplot(DT, aes(x=target, y=Ct, fill=sample)) + geom_bar(stat="identity", position="dodge")
            } else {
              DT <- pcr@data[,.(sample, target, Ct_mean)]
              p <- ggplot(DT, aes(x=target, y=Ct_mean, fill=sample)) + geom_bar(stat="identity", position="dodge")
            }

            return(p)
          })


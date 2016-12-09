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
            function(pcr, raw=FALSE, minimal=FALSE, plottype=c("bar", "point"), reverseGroups=FALSE, ...) standardGeneric("plotRQ"))
#' @aliases plotRQ,plotRQ-method
#' @rdname plotRQ
setMethod("plotRQ",
          signature("qPCR"),
          function(pcr, ...) {

            plottype <- match.arg(plottype)

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
              DT$group <- factor(ifelse(grepl("ko", DT$sample, ignore.case=TRUE), "KO", "WT"))

              meanz <- cbind(rbind(data.table(aggregate(RQ ~ target, data=DT[group=="KO"], FUN=mean)),
                                   data.table(aggregate(RQ ~ target, data=DT[group=="WT"], FUN=mean))),
                             group=rep(c("KO", "WT"), each=nrow(DT)/length(levels(DT$sample))))

              if (reverseGroups == TRUE) {
                DT$group    <-  factor(DT$group, levels=rev(levels(DT$group)))
                meanz$group <- factor(meanz$group, levels=rev(levels(meanz$group)))
              }

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
                              legend.title = element_text(size=18, family="Helvetica"),
                              aspect.ratio = 0.5)

              if (length(pcr@design) > 0) {
                if (length(levels(pcr@design)) <= 3) {

                  # design is defined, color appropriately
                  palettes <- c("Blues", "Greens", "Purples")
                  design.table <- table(pcr@design)
                  color.values <- character()
                  for (i in 1:length(design.table)) color.values <- c(color.values, rev(brewer.pal(9, palettes[i]))[1:design.table[i]])
                  if (plottype=="point") {
                    p <- p + scale_color_manual(values=color.values)
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

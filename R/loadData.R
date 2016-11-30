# ---------------------------------------------------------------------------- #
#' readPCR
#'
#' read raw qPCR data from an XLSX file
#'
#' details
#'
#' @param filename
#'
#' @importFrom xlsx read.xlsx
#' @importFrom data.table data.table setnames set
#' @export
setGeneric(
  name="readPCR",
  def=function(filename) {
    standardGeneric("readPCR")
  }
)
setMethod("readPCR",
          signature("character"),
          definition=function(filename) {
            raw.df <- read.xlsx(filename, 1, header=FALSE, stringsAsFactors=FALSE)
            meta1.dt <- data.table(raw.df[1:6,1:2]) # first 6 rows
            meta2.dt <- data.table(raw.df[(nrow(raw.df)-4):nrow(raw.df),1:2]) # last 5 rows
            pcr.dt   <- data.table(raw.df[(which(raw.df[,1]=="Well")+1):(nrow(raw.df)-4),])
            # first 28 column names should be stable, if there is more just copy from raw
            column.names <- c("well", "sample", "target", "task", "reporter", "quencher",
                              "RQ", "RQmin", "RQmax",
                              "Ct", "Ct_mean", "Ct_stdev",
                              "dCt", "dCt_mean", "dCt_SE", "HK_ctrl_dCt_mean", "HK_ctrl_dCt_SE", "ddCt",
                              "auto_Ct_thr", "Ct_thr",
                              "auto_baseline", "baseline_start", "baseline_end", "efficiency", "Tm1", "Tm2", "Tm3", "comments")
            if (ncol(raw.df) > 28) {
              column.names <- c(column.names, as.character(unname(unlist(raw.df[which(raw.df[,1]=="Well"), 29:ncol(raw.df)]))))
            }
            setnames(pcr.dt, column.names)
            # set proper data types for factor and numeric columns
            for (col in c(7L:15L, 18L, 20L, 25L, 26L, 27L)) set(pcr.dt, j=col, value=as.numeric(pcr.dt[[col]]))
            for (col in 1L:3L) set(pcr.dt, j=col, value=as.factor(pcr.dt[[col]]))

            pcr.list <- new("qPCR", data     = pcr.dt,
                                    metadata = rbind(meta1.dt, meta2.dt),
                                    pdata    = data.table())
            #pcr.list <- list(data     = pcr.dt,
            #                 metadata = rbind(meta1.dt, meta2.dt)) # for now, keep the data in a list with two data.tables (meta and data)

            return(pcr.list)
          }
)

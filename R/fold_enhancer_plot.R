fold_enhancer_plot.DESeqResults <- function(object, alpha = 0.05, bw = 0.01, breaks = waiver(), foldchange =0.5, yadjust=100, ladjust=0.25, radjust=1, xmin, xmax, ymax, main = deparse(substitute(object))) {
  if (is.null(object$log2FoldChange)) {
    stop("log2FoldChange column is not present: you should first run results()")
  }

  if (missing(xmin)) {
    xmin <- NA
  }

  if (missing(xmax)) {
    xmax <- max(object$log2FoldChange)
  }

  if (missing(ymax)) {
    ymax = 0.1*nrow(object)
  }

  FoldChange <- padj <- NULL
  isSig <- object[["padj"]] <= alpha
  df <- data.frame(padj = object[["padj"]],
                   FoldChange = 2^ object[["log2FoldChange"]],
                   isSig = isSig)

  ggplot() +
    geom_histogram(data = df, mapping = aes(x=FoldChange, fill = "grey"), alpha =0.5, color = "grey", binwidth =bw) +
    geom_histogram(data = subset(df, padj<=0.05 & FoldChange >= 1+ foldchange), mapping = aes(x=FoldChange, fill = "#2b8cbe"), alpha = 0.5, color = "#2b8cbe", binwidth =bw) +
    geom_histogram(data = subset(df, padj<=0.05 & (FoldChange < 1+ foldchange & FoldChange>=1)), mapping = aes(x=FoldChange, fill = "#a6bddb"), alpha = 0.5, color = "#a6bddb", binwidth =bw) +
    theme_classic() + scale_fill_identity("", guide = 'legend',labels = c(paste0(expression(p[FDR]<0.05)," & ", "enVars" ), paste0(expression(p[FDR]<0.05)," & ", "Non-enVars" ), "All Variants")) + xlab("Fold Change") + ylab("Count") +
    geom_vline(xintercept = 1+ foldchange, alpha = 0.4, color = "gold", size = 1)+ geom_text(aes(x=2+foldchange+radjust, label=paste(foldchange*100, "%\nchange"), y=yadjust), colour="black", size=4.5) +
    theme(axis.text=element_text(size=10),axis.title = element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10), legend.text=element_text(size=10), legend.position = c(0.75,0.3)) +
    annotation_logticks(base = 2, sides = "b") + annotation_logticks(base = 2, sides = "l") +
    scale_x_continuous(trans = "log2", limits = c(xmin,xmax), breaks = breaks) +
    scale_y_continuous(trans = "log2", limits = c(NA,ymax)) + ggtitle(main)
}

#' Plot histogram for significant allele with all alleles as background
#'
#' A quick plot to visualize the distribution of the significant alleles with all alleles as background. Yellow lines are added to indicate the fold change cut off.
#' @docType methods
#' @rdname fold_enhancer_plot
#' @aliases fold_enhancer_plot.DESeqResults
#' @param object a \code{DESeqResults} object produced by \code{\link{results}}
#' @param alpha the significance level for thresholding adjusted p-values
#' @param bw bandwidth
#' @param breaks One of:
#'   - `NULL` for no breaks
#'   - `waiver()` for the default breaks computed by the transformation object
#'   - A numeric vector of positions
#'   - A function that takes the limits as input and returns breaks as output
#' @param foldchange the fold change of the yellow line
#' @param yadjust the y location of fold change text
#' @param ladjust the adjustment for the left yellow line
#' @param radjust the adjustment for the right yellow line
#' @param xmin the minimum value for the x-axis
#' @param xmax the maximum value for the x-axis
#' @param ymax the maximum value for the y-axis
#' @param main the title of the figure.
#' @param ... additional arguments
#'
#' @importFrom ggplot2  ggplot aes annotation_logticks element_text geom_histogram geom_text geom_vline ggplot ggtitle scale_fill_identity scale_x_continuous scale_y_continuous theme theme_classic xlab ylab waiver
#' @export
#' @examples
#' dds <- DESeqDataSetFromMatrix()
#' dds <- DESeq(dds)
#' res <- results(dds, contrast=c("condition","exp","ctrl"))
#' fold_enhancer_plot(res)
setMethod("fold_enhancer_plot", signature(object="DESeqResults"), fold_enhancer_plot.DESeqResults)
setMethod("fold_enhancer_plot", signature(object="data.frame"), fold_enhancer_plot.DESeqResults)

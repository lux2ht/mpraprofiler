#' Dot plot for significant allelic different variant
#' @docType methods
#' @rdname allelic_enhancer_dot_plot
#' @param data data used for plotting
#' @param log2 column indicating the log2 fold change, character
#' @param max_fold column indicating maximum fold change, character
#' @param label column indicating the significant variant, character
#' @param foldchange fold change cut-off for claiming allelic enhancer variant, such as 0.25, numeric
#' @param type whether to include the annotation line with 1 means yes
#' @param main the tile of the graph
#'
#' @importFrom ggplot2  ggplot geom_point scale_color_manual scale_alpha_manual theme_classic scale_x_continuous scale_y_continuous geom_vline geom_hline geom_text xlab ylab theme ggtitle  guide_legend element_blank
#'
#' @export
#'
#' @examples
#' allelic_enhancer_dot_plot(te7_ctrl_plot,"log2_aver", "max_fold", "pos")
allelic_enhancer_dot_plot_1.25 <- function(data, log2, max_fold, label, foldchange = 0.25, type = 1, main = deparse(substitute(data))) {
  nr_r_ratio <- max_fold_change <- pos <- NULL
  data_plot <- as.data.frame(data)
  data_plot$nr_r_ratio <- 2^data_plot[, log2]
  data_plot$max_fold_change <- 2^data_plot[, max_fold]
  data_plot$pos <- as.factor(data_plot[, label])
  if (type==1) {
    ggplot() +
      geom_point(data = data_plot, mapping = aes(x=max_fold_change,y=nr_r_ratio, color = pos, alpha = pos)) +
      scale_color_manual(values=c("#999999", '#fb6a4a',"#de2d26"), guide = guide_legend(reverse=TRUE, override.aes = list(size = 2)),labels = c("Non-significant enVars", "Significant, non-allelic enVars", "Allelic enVars")) +
      scale_alpha_manual(values = c(0.1,1, 1), guide = FALSE) +
      #expression(p[FDR]<0.05)
      theme_classic()+
      scale_x_continuous(trans = "log2", breaks = c(0.5,1,2,4,8,16)) +
      scale_y_continuous(trans = "log2") +
      geom_vline(xintercept = 1,  color = "black") +
      geom_hline(yintercept = 2^log2(1+foldchange),  color = "gold", size = 1) +
      geom_hline(yintercept = 2^-log2(1+foldchange),  color = "gold", size = 1) +
      geom_hline(yintercept = 1,  color = "black") +
      geom_text(aes(x=2^-log2(2), label=paste0(foldchange*100, "% change"), y=2^log2(1+foldchange + 0.3)), colour="black", size=4) +
      geom_text(aes(x=2^-log2(2), label=paste( "-",foldchange*100, "% change"), y=2^-log2(1+foldchange + 0.3)), colour="black", size=4) +
      xlab("Enhancer fold change") + ylab("Genotype-dependence \n (Normalized non-reference\\reference fold change)") + ggtitle(main) +
      theme(axis.text=element_text(size=10),axis.title = element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_blank(), legend.text = element_text(size = 10), legend.position = c(0.8,0.8))
  } else {
    ggplot() +
      geom_point(data = data_plot, mapping = aes(x=max_fold_change,y=nr_r_ratio, color = pos, alpha = pos)) +
      scale_color_manual(values=c("#999999", '#fb6a4a',"#de2d26"), guide = guide_legend(reverse=TRUE, override.aes = list(size = 2)),labels = c("Non-significant enVars", "Significant, non-allelic enVars", "Allelic enVars")) +
      scale_alpha_manual(values = c(0.1,1, 1), guide = FALSE) +
      #expression(p[FDR]<0.05)
      theme_classic()+
      scale_x_continuous(trans = "log2", breaks = c(0.5,1,2,4,8,16)) +
      scale_y_continuous(trans = "log2") +
      geom_vline(xintercept = 1,  color = "black") +
      geom_hline(yintercept = 1,  color = "black") +
      ggtitle(main) +
      theme(axis.text=element_text(size=10),axis.title = element_text(size=10), plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_blank(), legend.text = element_text(size = 10), legend.position = c(0.8,0.8))
  }

}

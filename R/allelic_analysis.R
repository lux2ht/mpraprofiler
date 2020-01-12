#' A quickway to get the p value for each variant, only suited for one sample t-test
#'
#' @rdname p_value_function
#' @param x a set of data for calculating p-value
#'
#' @return p-value
#' @export
#'
p_value_function <- function(x) {
  p_value <- stats::t.test(x)$p.value
  return(p_value)
}


#' Get the max fold change for each varaint
#'
#' @rdname max_fold
#' @param object a \code{DESeqResults} object produced by \code{\link{results}}
#' @import dplyr
#' @return a data.frame containing variant name, max fold change and max absolute fold change
#' @export
#'
max_fold <- function(object){
  reference_id <- log2FoldChange <- snp <- max_fold_abs <- NULL
  #deal with the oligo_ID
  as.data.frame(object) %>% tibble::rownames_to_column("reference_id") %>% tidyr::separate(reference_id, c("snp", "nr/r", "allele"), sep = "_") -> object
  #pos
  object %>% dplyr::filter(log2FoldChange >=0) %>% dplyr::group_by(snp) %>% dplyr::summarise(max_fold = base::max(log2FoldChange), max_fold_abs = base::max(base::abs(log2FoldChange))) -> object_fold_pos
  #neg
  object %>% dplyr::filter(log2FoldChange <0) %>% dplyr::group_by(snp) %>% dplyr::summarise(max_fold = base::min(log2FoldChange), max_fold_abs = base::max(base::abs(log2FoldChange))) -> object_fold_neg
  #combined
  object_fold <- base::rbind(object_fold_pos, object_fold_neg)

  object_fold %>% dplyr::group_by(snp) %>% dplyr::summarise(max_fold=max_fold[which.max(max_fold_abs)]) -> object_fold_sum
  return(object_fold_sum)
}

#' Get the non-ref vs ref ratio for each variant
#'
#' @rdname nr_r_ratio
#' @param cts data counts prepared for DESeq2 analysis.
#' @import dplyr
#' @return a data.frame containing all the non-ref vs ref ratio with the correpsonding annotation
#' @export
#'
nr_r_ratio <- function(cts) {
  `reference_id` <- `nr/r` <- `snp` <- `.`<- NULL
  as.data.frame(cts) %>% tibble::rownames_to_column("reference_id") %>% tidyr::separate(reference_id, c("snp", "nr/r", "allele"), sep = "_") -> cts_sep

  cts_sep %>% dplyr::filter(`nr/r`=="Ref") -> cts_ref
  cts_sep %>% dplyr::filter(grepl("Non-Ref",`nr/r`)) -> cts_non_ref

  cts_join <- dplyr::inner_join(cts_ref,cts_non_ref, by= "snp", suffix = c(".ref", ".nonref"))

  cts_data <- cts_join %>% dplyr::mutate(snp2 = snp) %>% tidyr::unite("snp2","allele.nonref","allele.ref","nr/r.nonref",col = "compare_nrvsr", sep = "_") %>% dplyr::select(-"nr/r.ref")

  cts_nr_r_ratio <- (cts_data[,base::grepl("\\.nonref", base::colnames(cts_data))]/cts_data[,base::grepl("\\.ref", base::colnames(cts_data))]) %>% dplyr::bind_cols(cts_data[,c("snp", "compare_nrvsr")]) %>% stats::setNames(base::gsub(".nonref", "_nr/r_ratio", names(.)))

  return(cts_nr_r_ratio)
}

#' Get the experiment vs control comparision for all variants
#'
#' @param exp the experimental non-ref vs ref ratio, data.frame
#' @param ctrl the control non-ref vs ref ratio, data.frame
#' @param anno the annotation for each comparision, data.frame
#' @import dplyr
#' @return a data.frame containing all xperiment vs control comparision with the correpsonding annotation
#' @export
#'
allelic_compare <- function(exp, ctrl, anno) {
  data_norm <- exp/ctrl
  data_log2 <- data_norm  %>% dplyr::transmute_all(list(log2 = log2))
  data_log2$log2_aver <- base::apply(data_log2, 1, mean)
  data_combine <- base::cbind(data_norm,data_log2, anno)
  return(data_combine)
}


#' Get the allelic comparison p value for the enhancers
#'
#' @param allelic_compare the results produced by \code{\link{allelic_compare}}
#' @param enhancer a list of hte enhancer variant ID
#' @import dplyr
#' @return a data.frame containing all the allelic comparison p value and p.adjust for the enhancers
#' @export
#'
allelic_compare_p <- function(allelic_compare, enhancer) {
  snp <- NULL
  allelic_p <- allelic_compare %>% dplyr::filter(snp %in% enhancer)
  allelic_p$p_value <- base::apply(allelic_p[,base::grepl("_log2", names(allelic_p))],1,p_value_function)
  allelic_p$pFDR <- stats::p.adjust(allelic_p$p_value, method = "fdr")
  return(allelic_p)
}

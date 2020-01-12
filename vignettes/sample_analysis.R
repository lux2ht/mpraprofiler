# install.packages("devtools")
# install.packages("getPass")
devtools::install_github("lux2ht/mpraprofiler", force = T)
library(mpraprofiler)
library(readxl)
library(tidyverse)
library("DESeq2")
###############################prep the data###############################
#read all the data
file_list <- list.files("./data", pattern = '*.xlsx')
file_loc <- paste0("./data/", file_list)
df_list <- lapply(file_loc, read_xlsx)

#give each list a name
mpra_name <- str_split(file_list, "_", simplify = TRUE)[2:11,1]
names(df_list) <- c("EoE_ctrl", mpra_name)

#list to data frame
df <- bind_rows(df_list, .id = "id")

#select the useful ones for later analysis
df_use <- df %>% select(1,3,8)
names(df_use) <- c("sample","reference_id","counts")

#change to count matrix
count_table <- spread(df_use,sample, counts)
cts <- count_table

###############################Deseq2###############################
#build deseq2 matrix
cts <- cts[,2:ncol(cts)]
row.names(cts) <- count_table$reference_id
cts <- as.matrix(cts)

#design matrix
coldata <- data.frame("condition" = c("eoe_control", "te7","te7","te7","te7","te7","te7_il13","te7_il13","te7_il13","te7_il13","te7_il13"), row.names = colnames(cts))
#check the design matrix has the same order as the count matrix.
all(rownames(coldata) == colnames(cts))

#deseq2 analysis
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)

#te7 vs ctrl
te7_ctrl <- results(dds, contrast=c("condition","te7","eoe_control"))
#te7 with il13 stimulation vs ctrl
te7_il13 <- results(dds, contrast=c("condition","te7_il13","eoe_control"))
#te7 with il13 stimulation vs te7
te7_il13_te7 <- results(dds, contrast=c("condition","te7_il13","te7"))

###############################Enhancer###############################
#summary output
sink("te7_vs_ctrl_summary.txt")
summary(results(dds, contrast=c("condition","te7","eoe_control"),alpha = 0.05,lfcThreshold = log2(1.5)))
sink()

sink("te7_il13_vs_ctrl_summary.txt")
summary(results(dds, contrast=c("condition","te7_il13","eoe_control"),alpha = 0.05,lfcThreshold = log2(1.5)))
sink()

sink("te7_il13_vs_te7_summary.txt")
summary(results(dds, contrast=c("condition","te7_il13","te7"),alpha = 0.05,lfcThreshold = log2(1.15)))
sink()

#output the deseq2 result
write.csv(as.data.frame(te7_ctrl), "te7_vs_ctrl_deseq2.csv")
write.csv(as.data.frame(te7_il13), "te7_il13_vs_ctrl_deseq2.csv")
write.csv(as.data.frame(te7_il13_te7), "te7_il13_vs_te7_deseq2.csv")

###############################Plot###############################
#source the fold_enhancer_plot.R
fold_enhancer_plot(te7_ctrl, main = "TE7 vs Ctrl")
fold_enhancer_plot(te7_il13, main = "TE7 with IL13 stimulation vs Ctrl")
fold_enhancer_plot(te7_il13_te7, main = "TE7 with IL13 stimulation vs TE7",  xmin = 2^-log2(1.7), xmax = 2, bw = 0.001, breaks = c(1,2), ladjust = 0, radjust = 0, foldchange = 0.15)

###############################Histone and Homer###############################
#te7 enhancer
te7_ctrl_res <- as.data.frame(te7_ctrl)
te7_ctrl_res$reference_id <- rownames(te7_ctrl_res)
te7_ctrl_enhancer <- te7_ctrl_res %>% filter(log2FoldChange>=log2(1.5)&padj<=0.05)
te7_ctrl_snp <- unique(str_split(te7_ctrl_enhancer$reference_id, "_", simplify = TRUE)[,1])
write_lines(te7_ctrl_snp, "te7_ctrl_enhancer.txt", na = "NA", append = FALSE)

te7_ctrl_neg <- te7_ctrl_res %>% filter((log2FoldChange<=log2(1.1)&log2FoldChange>=-log2(1.1))&padj>0.05)
te7_ctrl_neg_snp <- unique(str_split(te7_ctrl_neg$reference_id, "_", simplify = TRUE)[,1])
write_lines(te7_ctrl_neg_snp, "te7_ctrl_neg.txt", na = "NA", append = FALSE)


#te7 with IL13 enhancer
te7_il13_res <- as.data.frame(te7_il13)
te7_il13_res$reference_id <- rownames(te7_il13_res)
te7_il13_enhancer <- te7_il13_res %>% filter(log2FoldChange>=log2(1.5)&padj<=0.05)
te7_il13_snp <- unique(str_split(te7_il13_enhancer$reference_id, "_", simplify = TRUE)[,1])
write_lines(te7_il13_snp, "te7_il13_enhancer.txt", na = "NA", append = FALSE)

te7_il13_neg <- te7_il13_res %>% filter((log2FoldChange<=log2(1.1)&log2FoldChange>=-log2(1.1))&padj>0.05)
te7_il13_neg_snp <- unique(str_split(te7_il13_neg$reference_id, "_", simplify = TRUE)[,1])
write_lines(te7_il13_neg_snp, "te7_il13_neg.txt", na = "NA", append = FALSE)

#te7 with IL13 enhancer vs te7
te7_il13_te7_res <- as.data.frame(te7_il13_te7)
te7_il13_te7_res$reference_id <- rownames(te7_il13_te7_res)

#either enhancers
either_enhancer <- full_join(te7_ctrl_enhancer,te7_il13_enhancer, by="reference_id")

#IL13 enhancer dependent enhancers p<=0.05 and 15%change
il13_enhancer <- te7_il13_te7_res %>% filter((reference_id %in% either_enhancer$reference_id) & padj<=0.05 & log2FoldChange>=log2(1.15))

il13_enhancer_snp <- unique(str_split(il13_enhancer$reference_id, "_", simplify = TRUE)[,1])
write_lines(il13_enhancer_snp, "il13_enhancer.txt", na = "NA", append = FALSE)

il13_enhancer_neg <- te7_il13_te7_res %>% filter((log2FoldChange<=log2(1.05)&log2FoldChange>=-log2(1.05))&padj>0.05)
il13_enhancer_neg_snp <- unique(str_split(il13_enhancer_neg$reference_id, "_", simplify = TRUE)[,1])
write_lines(il13_enhancer_neg_snp, "il13_enhancer_neg.txt", na = "NA", append = FALSE)

###############################Histone and Homer analysis###############################
#use ez_pipeline for analysis.


###############################Correlation###############################
normalized_cts <- counts(dds, normalized=TRUE)
cor_cts <- cor(normalized_cts)


###############################allelic analysis###############################
#add 0.5 to avoid infinite
cts <- cts +0.5
cts_nr_r_ratio <- nr_r_ratio(cts)
###############################Analysis###############################
###############################te7_vs_ctrl###############################
#all comparison 
#allelic dependence
#normalization and log2 transformation
te7_allelic_all <- allelic_compare(cts_nr_r_ratio[,2:6], cts_nr_r_ratio[,1], cts_nr_r_ratio[,c("snp","compare_nrvsr")])

#max fold change for two alleles
te7_ctrl_max_fold <- max_fold(te7_ctrl)

#do the p-value for the significant enhancer
te7_allelic_p <- allelic_compare_p(te7_allelic_all, te7_ctrl_snp)

write.csv(as.data.frame(te7_allelic_p), "te7_vs_ctrl_enhancer_allelic.csv")

#sig allielic snp
te7_allelic <- te7_allelic_p %>% filter(pFDR <= 0.05 & (log2_aver >= log2(1.20) | log2_aver <= -log2(1.20)))

write.csv(as.data.frame(te7_allelic), "te7_vs_ctrl_allelic.csv")

te7_allelic_snp <- te7_allelic %>% pull(snp) %>% unique()
write_lines(te7_allelic_snp, "te7_vs_ctrl_allelic_snp.txt", na = "NA", append = FALSE)

te7_allelic_snp_neg <- te7_allelic_p %>% filter(pFDR > 0.05 & (log2_aver < log2(1.05) & log2_aver > -log2(1.05))) %>% pull(snp) %>% unique()
write_lines(te7_allelic_snp_neg, "te7_vs_ctrl_allelic_neg.txt", na = "NA", append = FALSE)

#combined together for plot
te7_ctrl_plot <- inner_join(te7_allelic_all,te7_ctrl_max_fold, by = "snp") %>% left_join(te7_allelic_p[,c("snp","p_value", "pFDR")], by = "snp")

te7_ctrl_plot$pos <- ifelse(te7_ctrl_plot$snp %in% te7_allelic_snp & (te7_ctrl_plot$log2_aver >= log2(1.20) | te7_ctrl_plot$log2_aver <= -log2(1.20)), 1,0) 
te7_ctrl_plot %>% replace_na(list(pos=0)) ->te7_ctrl_plot

allelic_enhancer_dot_plot(te7_ctrl_plot, log2 = "log2_aver", max_fold = "max_fold", label = "pos")
allelic_enhancer_dot_plot(te7_ctrl_plot, log2 = "log2_aver", max_fold = "max_fold", label = "pos", type = 0)

###############################te7_il13_vs_ctrl###############################
#all comparison 
#allelic dependence
#normalization and log2 transformation
te7_il13_allelic_all <- allelic_compare(cts_nr_r_ratio[,7:11],cts_nr_r_ratio[,1],(cts_nr_r_ratio[,c("snp","compare_nrvsr")]))

#max fold change for two alleles
te7_il13_max_fold <- max_fold(te7_il13)

#do the p-value for the significant enhancer
te7_il13_allelic_p <- allelic_compare_p(te7_il13_allelic_all,te7_il13_snp)

#sig allielic snp
te7_il13_allelic <- te7_il13_allelic_p %>% filter(pFDR <= 0.05 & (log2_aver >= log2(1.20) | log2_aver <= -log2(1.20)))

write.csv(as.data.frame(te7_il13_allelic), "te7_il13_vs_ctrl_allelic.csv")

te7_il13_allelic_snp <- te7_il13_allelic %>% pull(snp) %>% unique()
write_lines(te7_il13_allelic_snp, "te7_il13_vs_ctrl_allelic_snp.txt", na = "NA", append = FALSE)

te7_il13_allelic_snp_neg <- te7_il13_allelic_p %>% filter(pFDR > 0.05 & (log2_aver < log2(1.05) & log2_aver > -log2(1.05))) %>% pull(snp) %>% unique()
write_lines(te7_il13_allelic_snp_neg, "te7_il13_vs_ctrl_allelic_neg.txt", na = "NA", append = FALSE)

#combined together for plot
te7_il13_plot <- inner_join(te7_il13_allelic_all,te7_il13_max_fold, by = "snp") %>% left_join(te7_il13_allelic_p[,c("snp","p_value", "pFDR")], by = "snp")

te7_il13_plot$pos <- ifelse(te7_il13_plot$pFDR <=0.05, 1,0) 
te7_il13_plot %>% replace_na(list(pos=0)) ->te7_il13_plot

allelic_enhancer_dot_plot(te7_il13_plot, log2 = "log2_aver", max_fold = "max_fold", label = "pos")



###############################te7_il13_vs_te7###############################
#all comparison 
#allelic dependence
#normalization and log2 transformation
te7_il13_te7_allelic_all <- allelic_compare(cts_nr_r_ratio[,7:11],cts_nr_r_ratio[,2:6],cts_nr_r_ratio[,c("snp","compare_nrvsr")])

#max fold change for two alleles
te7_il13_te7_max_fold <- max_fold(te7_il13_te7)

#do the p-value for the significant enhancer
te7_il13_te7_allelic_p <- allelic_compare_p(te7_il13_te7_allelic_all,il13_enhancer_snp)

#sig allielic snp ## No fold change select
te7_il13_te7_allelic <- te7_il13_te7_allelic_p %>% filter(pFDR <= 0.05)

write.csv(as.data.frame(te7_il13_te7_allelic), "te7_il13_vs_te7_allelic.csv")

te7_il13_te7_allelic_snp <- te7_il13_te7_allelic %>% pull(snp) %>% unique()
write_lines(te7_il13_te7_allelic_snp, "te7_il13_vs_te7_allelic_snp.txt", na = "NA", append = FALSE)
#rs73131258
#https://eqtl.onderzoek.io/index.php?page=gene_cis_details&gene=RNF114
#https://www.ncbi.nlm.nih.gov/pubmed/28165122
#rs6875763
#http://eqtl.rc.fas.harvard.edu/eqtlbrowser/mrcau133list/42267
#http://europepmc.org/abstract/MED/23345460

te7_il13_te7_allelic_snp_neg <- te7_il13_te7_allelic_p %>% filter(pFDR > 0.05) %>% pull(snp) %>% unique()
write_lines(te7_il13_te7_allelic_snp_neg, "te7_il13_vs_te7_allelic_neg.txt", na = "NA", append = FALSE)

#combined together for plot
te7_il13_te7_plot <- inner_join(te7_il13_te7_allelic_all,te7_il13_te7_max_fold, by = "snp") %>% left_join(te7_il13_te7_allelic_p[,c("snp","p_value", "pFDR")], by = "snp")

te7_il13_te7_plot$pos <- ifelse(te7_il13_te7_plot$pFDR <=0.05, 1,0) 
te7_il13_te7_plot %>% replace_na(list(pos=0)) ->te7_il13_te7_plot

allelic_enhancer_dot_plot(te7_il13_te7_plot, log2 = "log2_aver", max_fold = "max_fold", label = "pos", type = 0)

###############################summary###############################
snp_annation <- as.data.frame(unique(cts_nr_r_ratio$snp), stringsAsFactors=FALSE)
colnames(snp_annation) <- "SNP_ID"
snp_annation$Enhancer_activity_without_stimulation <- ifelse(snp_annation$SNP_ID %in% te7_ctrl_snp, 1, 0)
snp_annation$Enhancer_activity_with_stimulation <- ifelse(snp_annation$SNP_ID %in% te7_il13_snp, 1, 0)
snp_annation$Allelic_activity_without_stimulation <- ifelse(snp_annation$SNP_ID %in% te7_allelic_snp, 1, 0)
snp_annation$Allelic_activity_with_stimulation <- ifelse(snp_annation$SNP_ID %in% te7_il13_allelic_snp, 1, 0)
snp_annation$Stimulation_affected_enhancer_activity  <- ifelse(snp_annation$SNP_ID %in% il13_enhancer_snp, 1, 0)
snp_annation$Stimulation_affected_allelic_activity  <- ifelse(snp_annation$SNP_ID %in% te7_il13_te7_allelic_snp, 1, 0)

write.csv(snp_annation, "snp_annotation.csv")

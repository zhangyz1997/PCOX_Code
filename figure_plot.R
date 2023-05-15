setwd(this.path::this.dir())
## PCOX 正式出图代码
### 先准备数据
if (T) {
  library("this.path")
  library(IRanges)
  library(IOBR)
  library(tidyverse)
  setwd(this.dir()) # https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
  library(readxl)
  library(writexl)
  library(maftools)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggpubr)
  library(msigdbr)
  library(ReactomePA)
  library(tinyarray)
  library(patchwork)
  library(cowplot)
  excel2maf <- function(wes_result) {
    wes_maf <- data.frame(
      Hugo_Symbol = wes_result$symbol,
      Chromosome = stringr::str_sub(wes_result$chr, 4),
      Start_Position = as.numeric(wes_result$position),
      End_Position = as.numeric(wes_result$position),
      Tumor_Sample_Barcode = wes_result$patient_id,
      Reference_Allele = wes_result$ref,
      Tumor_Seq_Allele2 = wes_result$alt,
      # Variant_Type = "SNP",
      HGVSp_Short = wes_result$aminochange_ts,
      VAF = wes_result$tumor_alt_depth / (
        wes_result$tumor_alt_depth + wes_result$tumor_refer_depth
      ) * 100
    ) %>%
      mutate(
        # 参考 https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl 转化突变类型
        Variant_Classification = case_when(
          wes_result$mutationtype %in% c(
            "splice_acceptor_variant",
            "splice_donor_variant",
            "transcript_ablation",
            "exon_loss_variant"
          ) ~ "Splice_Site",
          wes_result$mutationtype == "splice_region_variant" ~ "Splice_Region",
          wes_result$mutationtype == "3_prime_UTR_variant" ~ "3'UTR",
          wes_result$mutationtype == "5_prime_UTR_variant" ~ "5'UTR",
          wes_result$mutationtype == "downstream_gene_variant" ~ "3'Flank",
          wes_result$mutationtype == "upstream_gene_variant" ~ "5'Flank",
          wes_result$mutationtype == "intron_variant" ~ "Intron",
          wes_result$mutationtype == "missense_variant" ~ "Missense_Mutation",
          wes_result$mutationtype == "non_coding_transcript_exon_variant" ~ "RNA",
          wes_result$mutationtype == "start_lost" ~ "Translation_Start_Site",
          wes_result$mutationtype == "stop_lost" ~ "Nonstop_Mutation",
          wes_result$mutationtype == "stop_gained" ~ "Nonsense_Mutation",
          wes_result$mutationtype %in% c("stop_retained_variant", "synonymous_variant") ~ "Silent",
          wes_result$mutationtype == "frameshift_variant" &
            stringr::str_detect(wes_result$cdnachange_ts, "del") ~ "Frame_Shift_Del",
          wes_result$mutationtype == "frameshift_variant" &
            (
              stringr::str_detect(wes_result$cdnachange_ts, "ins") |
                stringr::str_detect(wes_result$cdnachange_ts, "dup")
            ) ~ "Frame_Shift_Ins",
          wes_result$mutationtype == "inframe_deletion" ~ "In_Frame_Del",
          wes_result$mutationtype == "inframe_insertion" ~ "In_Frame_Ins"
        ),
        Variant_Type = case_when(
          wes_result$mutationtype == "frameshift_variant" &
            stringr::str_detect(wes_result$cdnachange_ts, "del") ~ "DEL",
          wes_result$mutationtype == "frameshift_variant" &
            (
              stringr::str_detect(wes_result$cdnachange_ts, "ins") |
                stringr::str_detect(wes_result$cdnachange_ts, "dup")
            ) ~ "INS",
          wes_result$mutationtype == "inframe_deletion" ~ "DEL",
          wes_result$mutationtype == "inframe_insertion" ~ "INS",
          TRUE ~ "SNP"
        )
      )
    return(wes_maf)
  }
  clinic <- read_xlsx("./data/clinical.xlsx") %>% filter(is.na(`编号`) == FALSE)
  clinic$Tumor_Sample_Barcode <- clinic$'编号'
  clinic$response[clinic$response == "/"] <- 1
  clinic$PD[clinic$PD == "/"] <- 0
  clinic <- clinic %>% filter(姓名 != "柯国浩")
  
  clinic_detailed <- read_xlsx("./data/20220524更新1005.xlsx") %>% filter(姓名 %in% clinic$姓名)
  clinic$surgery <- clinic_detailed$surgery[match(clinic$姓名, clinic_detailed$姓名)]
  clinic$cr <- clinic_detailed$研究者_CR[match(clinic$姓名, clinic_detailed$姓名)]
  
  immport_genelist <-
    read.table(file = "./data/GeneList.txt",
               header = T,
               sep = "\t")
  ## 读取WES的数据并整理为MAF格式
  wes_snp <- read_xlsx("./data/wes.xlsx", sheet = "SNP")
  wes_cnv <- read_xlsx("./data/wes.xlsx", sheet = "CNV")
  wes_maf <- read.maf(maf = "./data/group_all.maf")@data
  ### 专门读取柯国浩的WES数据
  wes_snp_kgh <-
    read.table(
      "./data/A14525_柯国浩/A14525_N5251803D_T5253418D_P011452501_TNscope_somatic_pass_vep_snp_filtered.txt",
      header = TRUE,
      sep = "\t"
    )
  wes_indel_kgh <-
    read.table(
      "./data/A14525_柯国浩/A14525_N5251803D_T5253418D_P011452501_TNscope_somatic_pass_vep_indel_filtered.txt",
      header = TRUE,
      sep = "\t"
    )
  
  wes_maf$Tumor_Sample_Barcode <-
    stringr::str_sub(wes_maf$Tumor_Sample_Barcode, 1, 6)
  wes_maf_snp <- data.frame(
    Hugo_Symbol = wes_snp$symbol,
    Chromosome = stringr::str_sub(wes_snp$chr, 4),
    Start_Position = as.numeric(wes_snp$position),
    End_Position = as.numeric(wes_snp$position),
    Tumor_Sample_Barcode = wes_snp$patient_id,
    Reference_Allele = wes_snp$ref,
    Tumor_Seq_Allele2 = wes_snp$alt,
    Variant_Type = "SNP",
    HGVSp_Short = wes_snp$aminochange_ts,
    VAF = wes_snp$tumor_alt_depth / (wes_snp$tumor_alt_depth + wes_snp$tumor_refer_depth) * 100
  ) %>%
    mutate(
      # 参考 https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl 转化突变类型
      Variant_Classification = case_when(
        wes_snp$mutationtype %in% c(
          "splice_acceptor_variant",
          "splice_donor_variant",
          "transcript_ablation",
          "exon_loss_variant"
        ) ~ "Splice_Site",
        wes_snp$mutationtype == "splice_region_variant" ~ "Splice_Region",
        wes_snp$mutationtype == "3_prime_UTR_variant" ~ "3'UTR",
        wes_snp$mutationtype == "5_prime_UTR_variant" ~ "5'UTR",
        wes_snp$mutationtype == "downstream_gene_variant" ~ "3'Flank",
        wes_snp$mutationtype == "upstream_gene_variant" ~ "5'Flank",
        wes_snp$mutationtype == "intron_variant" ~ "Intron",
        wes_snp$mutationtype == "missense_variant" ~ "Missense_Mutation",
        wes_snp$mutationtype == "non_coding_transcript_exon_variant" ~ "RNA",
        wes_snp$mutationtype == "start_lost" ~ "Translation_Start_Site",
        wes_snp$mutationtype == "stop_lost" ~ "Nonstop_Mutation",
        wes_snp$mutationtype == "stop_gained" ~ "Nonsense_Mutation",
        wes_snp$mutationtype %in% c("stop_retained_variant", "synonymous_variant") ~ "Silent"
      )
    ) %>%
    filter(Tumor_Sample_Barcode %in% wes_maf$Tumor_Sample_Barcode == FALSE)
  wes_maf_snp_kgh <- excel2maf(wes_snp_kgh)
  wes_maf_indel_kgh <- excel2maf(wes_indel_kgh)
  
  cntable_wes <- wes_cnv %>%
    dplyr::select(symbol, patient_id, haploid) %>%
    mutate(cnv_status = case_when(haploid >= 3 ~ "Amp",
                                  haploid <= 1.2 ~ "Del",
                                  TRUE ~ NA_character_)) %>%
    filter(is.na(cnv_status) == FALSE) %>%
    dplyr::select(symbol, patient_id, cnv_status) %>%
    unique()
  # construct_maf <- read.maf(wes_maf_snp, clinicalData = clinic, cnTable = cntable_wes)
  # construct_maf <- read.maf(wes_maf, clinicalData = clinic, cnTable = cntable_wes)
  # construct_maf <- read.maf(wes_maf, clinicalData = clinic)
  construct_maf <-
    merge_mafs(list(
      read.maf(wes_maf),
      read.maf(wes_maf_snp),
      read.maf(wes_maf_snp_kgh),
      read.maf(wes_maf_indel_kgh)
    ))
  construct_maf <-
    read.maf(
      construct_maf@data,
      clinicalData = clinic,
      cnTable = cntable_wes,
      rmFlags = TRUE
    )
}
response_enrich <- clinicalEnrichment(
  maf = construct_maf,
  clinicalFeature = "response",
  annotationDat = NULL,
  minMut = 2,
  useCNV = FALSE,
  pathways = FALSE
)
PFS_enrich <- clinicalEnrichment(
  maf = construct_maf,
  clinicalFeature = "PD",
  annotationDat = NULL,
  # minMut = 5,
  useCNV = FALSE,
  pathways = FALSE
)
plotEnrichmentResults(response_enrich,
                      pVal = 0.05,
                      ORthr = 1,
                      featureLvls = NULL,
                      cols = NULL,
                      # annoFontSize = 0.8,
                      # geneFontSize = 0.8,
                      # legendFontSize = 0.8
)

sig_genes <- response_enrich$groupwise_comparision %>% filter(p_value <= 0.05 & OR >= 1) %>% dplyr::select(Hugo_Symbol) %>% unique()
lollipopPlot2(m1 = subsetMaf(construct_maf, clinQuery = "response == 1"), m2 = subsetMaf(construct_maf, clinQuery = "response == 0"), gene = sig_genes[1], AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "Responders", m2_name = "Nonresponders", m1_label = "all", m2_label = "all")
 
### TMB相关图
maf.mutload <- maftools::tmb(maf = construct_maf, captureSize = 35, logScale = TRUE)
canc_tmb <- merge(x = data.frame(maf.mutload), y = clinic %>% dplyr::select(c( "编号","response","PD")),
                  by.x = "Tumor_Sample_Barcode",
                  by.y = "编号")
canc_tmb$response[which(canc_tmb$response == "/")] <- "1"
canc_tmb$PD[which(canc_tmb$PD == "/")] <- "0"
canc_tmb$response <- as.factor(canc_tmb$response)
canc_tmb$PD <- as.factor(canc_tmb$PD)
canc_tmb$BEST_RESPONSE <- case_when(
  canc_tmb$response == "1" ~ "CR/PR",
  canc_tmb$PD == "1" ~ "PD",
  TRUE ~ "SD"
)
writexl::write_xlsx(canc_tmb, path = "./figures/TMB.xlsx")

#### ggplot2重制
tmb_plot <- ggplot(data = canc_tmb %>% arrange(total_perMB_log), mapping = aes(x = 1:nrow(canc_tmb), y = total_perMB, group = BEST_RESPONSE))+
  geom_point(aes(color=BEST_RESPONSE))+
  scale_y_log10()+
  geom_hline(yintercept = median(canc_tmb[,"total_perMB"], na.rm = TRUE), linetype = "dashed")+
  labs(title = "Mutation Burden", x = "Patients", y = "TMB/MB", colour = "Response", caption = paste0("Median: ", round(median(maf.mutload[,total_perMB], na.rm = TRUE),3), "/MB"))+
  theme_classic()+
  # theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.position = c(0.8, 0.33))
  # scale_y_continuous(breaks = pretty(range(canc_tmb[canc_tmb$total != 0,"total_perMB_log"])))
plot(tmb_plot)

# ggplot(canc_tmb,aes(x = response,y = total_perMB_log)) + 
#   geom_boxplot() + 
#   labs(x = "Response", y = "logTMB") +
#   stat_compare_means(method = "wilcox.test")
violin_response <- ggplot(canc_tmb,aes(x = response,y = total_perMB)) + 
  # geom_violin() + 
  # geom_boxplot(width=0.15)+
  geom_boxplot(aes(colour = response))+
  geom_point(aes(colour = response))+
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), comparisons = list(c("0","1")))+
  # ggthemes::theme_clean()+
  scale_y_log10()+
  scale_x_discrete(labels = c("SD/PD","CR/PR"))+
  scale_color_discrete(labels = c("SD/PD","CR/PR"))+
  labs(x = "Best Response", y = "TMB/MB") +
  stat_compare_means(method = "wilcox.test",size = 2)+theme_classic()+
  theme(legend.position = "none")
violin_pd <- ggplot(canc_tmb,aes(x = PD,y = total_perMB)) + 
  # geom_violin() + 
  # geom_boxplot(width=0.15)+
  geom_boxplot(aes(colour = PD))+
  geom_point(aes(colour = PD))+
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), comparisons = list(c("0","1")))+
  # ggthemes::theme_clean()+
  scale_y_log10()+
  scale_x_discrete(labels = c("Non-PD","PD"))+
  scale_color_discrete(labels = c("Non-PD","PD"))+
  labs(x = "Tumor Progression", y = "TMB/MB") +
  stat_compare_means(method = "wilcox.test",size = 2)+theme_classic()+ 
  theme(legend.position = "none")
wilcox.test(x = canc_tmb$total_perMB[which(canc_tmb$response == 1)], y = canc_tmb$total_perMB[which(canc_tmb$response == 0)])

violins <- violin_response + violin_pd
final_tmb <- tmb_plot / violins + plot_layout(widths = c(2,1))
plot(final_tmb)
ggsave(filename = "tmb.pdf", plot = final_tmb, device = "pdf", path = "./figures",width = 150, height = 260/2, units = "mm")

### OncoPrint
pdf("./figures/oncoplot_1.pdf", width = 5.5, height = 4)
oncoplot(
  maf = construct_maf,
  # top = 10,
  clinicalFeatures = c("response","PD","cr","surgery"),
  genes = c("KRAS","NRAS","BRAF", "JAK1","JAK2", "B2M", "MDM2","MDM4","CTNNB1"),
  drawColBar = TRUE,
  logColBar = FALSE,
  drawRowBar = FALSE,
  topBarData = data.frame(canc_tmb$Tumor_Sample_Barcode, TMB = canc_tmb$total_perMB),
  includeColBarCN = FALSE,
  # fontSize = 0.6,
  legendFontSize = 0.9,
  annotationFontSize = 1.0,
  anno_height = 3,
  # pathways = "auto",
  sortByAnnotation = TRUE,
  cBioPortal = TRUE,
  removeNonMutated = FALSE,
  annoBorderCol = "white",
  borderCol = "white",
  bgCol = "white",
  drawBox = TRUE
)
dev.off()

## KEGG通路富集
library(org.Hs.eg.db)
library(clusterProfiler)
kegg_SYMBOL_hsa <- function(genes){ 
  gene.df <- bitr(genes, fromType = "SYMBOL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df) 
  diff.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 0.99,
                        qvalueCutoff = 0.99
  )
  return(setReadable(diff.kk, OrgDb = org.Hs.eg.db,keyType = 'ENTREZID'))
}

mut_matrix <- mutCountMatrix(construct_maf)
geneList <- list()
geneList[["response"]] <- mut_matrix[, intersect(colnames(mut_matrix), clinic$编号[clinic$response == "1" & is.na(clinic$编号) == FALSE])] %>% rowSums(.)
geneList[["response"]] <- names(geneList[["response"]][which(geneList[["response"]] != 0)])
geneList[["noresponse"]] <- mut_matrix[, intersect(colnames(mut_matrix), clinic$编号[clinic$response == "0" & is.na(clinic$编号) == FALSE])] %>% rowSums(.)
geneList[["noresponse"]] <- names(geneList[["noresponse"]][which(geneList[["noresponse"]] != 0)])

kegglist <- list()
kegglist[["response"]] <- kegg_SYMBOL_hsa(genes = geneList[["response"]])
barplot(kegglist[["response"]])
kegg1 <- dotplot(kegglist[["response"]])+ggtitle("Responders")
kegglist[["noresponse"]] <- kegg_SYMBOL_hsa(genes = geneList[["noresponse"]])
barplot(kegglist[["noresponse"]])
kegg2 <- dotplot(kegglist[["noresponse"]])+ggtitle("Non-Responders")

kegg_combined <- kegg1 / kegg2
ggsave(filename = "kegg_mutation.pdf", plot = kegg_combined, path = "./figures",width = 150, height = 200, units = "mm")

writexl::write_xlsx(list(response = kegglist[["response"]]@result, noresponse = kegglist[["noresponse"]]@result), path = "./figures/kegg_mutation.xlsx")

### 转录组出图
rna_seq_data <- read_xlsx(path = "./data/rnaseq.xlsx")
group_list_complete <- merge(
  rna_seq_data %>% dplyr::select(sample_id, patient_id),
  clinic %>% dplyr::select(c("姓名", "编号","PFS_event","response","PD")),
  by.x = "patient_id",
  by.y = "编号"
) %>% unique() %>%
  magrittr::set_rownames(.$sample_id) %>%
  dplyr::select(-"sample_id") %>% unique() %>%
  rownames_to_column(var = "sample_id") %>%
  magrittr::set_rownames(.$sample_id) %>% 
  mutate(benefit = ifelse((response == "1" & PD == "0"),"1","0"))
raw_count_matrix <- rna_seq_data %>%
  dplyr::select(sample_id, gene, count) %>%
  filter(sample_id %in% group_list_complete$sample_id) %>%
  unique() %>%
  pivot_wider(names_from = sample_id, values_from = count) %>%
  column_to_rownames(var = "gene")

## 构建TPM矩阵
tpm_matrix <- rna_seq_data %>%
  dplyr::select(sample_id, gene, tpm) %>%
  filter(sample_id %in% group_list_complete$sample_id) %>%
  unique() %>%
  pivot_wider(names_from = sample_id, values_from = tpm) %>%
  column_to_rownames(var = "gene")

## 获得log2用的标准化矩阵
tpm_log2_matrix <- rna_seq_data %>%
  mutate(tpm_log2 = log2(tpm+1)) %>%
  dplyr::select(sample_id, gene, tpm_log2) %>%
  filter(sample_id %in% group_list_complete$sample_id) %>%
  unique() %>%
  pivot_wider(names_from = sample_id, values_from = tpm_log2) %>%
  column_to_rownames(var = "gene")

library(IOBR)
group <- group_list_complete[colnames(raw_count_matrix),"response"] ## response为1,因此展示的是response相对于non-response的结果
group[group == "/"] <- "1"
project_id = "PCOX_Response"
## IOBR用的数据是 log2(TPM+1) 标准化后的表达矩阵 (tpm_log2_matrix), 且列名为样品
pdata_group <- data.frame(ID = colnames(tpm_log2_matrix), Response = as.factor(group))
which(apply(tpm_log2_matrix, 1, var)==0)
## 首先推算内建的所有的signature
sig_res<-calculate_sig_score(
  pdata           = NULL,
  eset            = tpm_log2_matrix[which(apply(tpm_log2_matrix, 1, var)!=0), ], # 除去全体样本中表达量一致(为0)的基因,避免 https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var 的报错
  signature       = signature_collection,
  method          = "ssGSEA",
  mini_gene_count = 2)
sig_hallmark<-calculate_sig_score(
  pdata           = NULL,
  eset            = tpm_log2_matrix[which(apply(tpm_log2_matrix, 1, var)!=0), ],
  signature       = hallmark,
  method          = "ssGSEA",
  mini_gene_count = 2)
sig_kegg<-calculate_sig_score(
  pdata           = NULL,
  eset            = tpm_log2_matrix[which(apply(tpm_log2_matrix, 1, var)!=0), ],
  signature       = kegg,
  method          = "ssGSEA",
  mini_gene_count = 2)
## 接着批量推算TME和免疫相关signature与response的关系
iobr_cor_plot(pdata_group           = pdata_group,
              id1                   = "ID",
              feature_data          = sig_res,
              id2                   = "ID",
              target                = NULL,
              group                 = "Response",
              is_target_continuous  = FALSE,
              padj_cutoff           = 1,
              index                 = 1,
              category              = "signature",
              signature_group       = sig_group,
              ProjectID             = project_id,
              palette_box           = "paired1",
              palette_corplot       = "pheatmap",
              palette_heatmap       = 2,
              feature_limit         = 20,
              character_limit       = 30,
              show_heatmap_col_name = TRUE,
              show_col              = FALSE,
              show_plot             = TRUE,
              path                  = paste0("./output/",project_id,"/"))


## DEG热图
library(readxl)
library(pheatmap)
library(ComplexHeatmap)
library(scales)
DEG <- read_excel("./output/PCOX_Response/DEGsPCOX_Response_edger.xlsx")
if("FDR" %in% colnames(DEG)){
  threshold<-as.factor((DEG$logFC>0.5|DEG$logFC<(-0.5)) & DEG$PValue<0.05)
}else if("log2FoldChange" %in% colnames(DEG)){
  threshold<-as.factor((abs(DEG$log2FoldChange)>=0.5) & DEG$pvalue<0.05)
}
DEG_sig <- as.data.frame(DEG[which(threshold == TRUE),])
rownames(DEG_sig) <- DEG_sig[,1]

deg_heatmap_matrix <- tpm_log2_matrix[rownames(DEG_sig), ] %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = !contains("gene"),
    names_to = "patients",
    values_to = "tpm",
    values_drop_na = TRUE
  ) %>%
  left_join(
    x = ., y = group_list_complete, by = c("patients" = "sample_id")
  ) %>% mutate(response = as.factor(response), PFS_event = as.factor(PFS_event), PD = as.factor(PD))
deg_heatmap <- deg_heatmap_matrix %>%
  group_by(response) %>%
  heatmap(
    .column = patients,
    .row = gene,
    .value = tpm,
    scale = "row",
    palette_value = c("red", "white", "blue"),
    palette_grouping = list(c(rgb(223,141,143, maxColorValue = 255),rgb(126,162,237, maxColorValue = 255))),
    show_row_names = FALSE,
    show_column_names = FALSE
  ) %>%
  add_tile(PD, palette = c("#d4e967","#d4c9ed")) %>%
  add_tile(response, palette = c("#89debd","#ff9663"))
print(deg_heatmap)
deg_heatmap %>% save_pdf("./figures/Figure_d_Response_new.pdf")

## 柱状图
bargenes <- c("IFNG", "GZMA", "GZMB", "IDO1", "LAG3", "PTEN", "TAPBP","PSME2")
barchart_matrix <- tpm_matrix[bargenes, ] %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = !contains("gene"),
    names_to = "patients",
    values_to = "tpm",
    values_drop_na = TRUE
  ) %>%
  left_join(
    x = ., y = group_list_complete, by = c("patients" = "sample_id")
  ) %>% mutate(response = as.factor(response), PFS_event = as.factor(PFS_event), PD = as.factor(PD))
barcharts <- lapply(bargenes, function(x){
  barchart_single <- ggplot(data = barchart_matrix %>% filter(gene == x), aes(x = response, y = tpm, colour = response)) +
    geom_boxplot()+
    geom_point(aes(colour = response))+
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), comparisons = list(c("0","1")))+
    ggthemes::theme_clean()+
    ggtitle(x)+
    labs(y = "TPM", colour = "Response")+
    scale_x_discrete(labels = c("SD/PD","CR/PR"))+
    scale_color_discrete(labels = c("SD/PD","CR/PR"))+
    theme(axis.title.x = element_blank())
  return(barchart_single)
})
library(patchwork)
barcharts_combined <- (patchwork::wrap_plots(barcharts, ncol = 4, nrow = 2, guides = "collect")) & theme(legend.position = "bottom")
barcharts_combined
ggsave(filename = "./figures/Figure_e_Barchart.pdf", barcharts_combined, width = 16, height = 14, units = "cm")

## ssGSEA细胞热图
library(GSVA)
immune_gmt <- read.csv("./data/ssGSEA_immune.csv")
immune_gmt_list <- split(as.matrix(immune_gmt)[,1], immune_gmt[,2])
gsva_immune <- gsva(as.matrix(tpm_log2_matrix), immune_gmt_list, mx.diff=FALSE, verbose=FALSE, parallel.sz=2, method = "ssgsea")

DEG_immune <- list()
divide_group <- function(division_criteria){
  if(division_criteria == "response"){
    group <- group_list_complete[colnames(raw_count_matrix),"response"] ## response为1,因此展示的是response相对于non-response的结果
    group[group == "/"] <- "1"
    project_id = "PCOX_Response"
  }else if(division_criteria == "pd"){
    group <- group_list_complete[colnames(raw_count_matrix),"PD"] ## response为1,因此展示的是response相对于non-response的结果
    group[group == "/"] <- "0"
    project_id = "PCOX_PD"
  }else if(division_criteria == "benefit"){
    group <- ifelse((group_list_complete$response == "1" & group_list_complete$PD == "0"),"1","0")
    project_id <- "PCOX_Benefit"
  }
  return(group)
}

division_criteria <- "response"
group <- divide_group(division_criteria)
library(limma)
design = model.matrix(~group)
fit = lmFit(gsva_immune, design)
fit = eBayes(fit)
DEG_immune[[division_criteria]] = topTable(fit, coef = 2, number = Inf)
head(DEG_immune[[division_criteria]],20)
tinyarray::draw_heatmap(gsva_immune[head(rownames(DEG_immune[[division_criteria]]),20),],group,show_rownames = TRUE, split_column = TRUE,show_column_title = TRUE, color = (colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlBu"))))(100), color_an = palettes(category = "box", palette = "nrc", show_col = FALSE, show_message = FALSE))

library(ComplexHeatmap)
library(scales)
library(RColorBrewer)

ssGSEA_mat_scaled <- (gsva_immune %>% pheatmap:::scale_rows() %>% as.data.frame %>% dplyr::select(clinic$RNA_FILE))
rownames(ssGSEA_mat_scaled) <- gsub(" CIBERSORT","",gsub("\\_"," ",rownames(ssGSEA_mat_scaled)))

rg <- range(ssGSEA_mat_scaled)
rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
rg[2] = rg[2] + (rg[2] - rg[1])* 0.02

anno_multiple_boxplot = function(index) {
  nr = length(index)
  pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
  for(i in seq_along(index)) {
    grid.rect(y = nr-i+1, height = 1, default.units = "native")
    grid.boxplot(unlist(ssGSEA_mat_scaled[index[i], which(group == 0)]), pos = nr-i+1 + 0.2, box_width = 0.3, 
                 gp = gpar(fill = 3), direction = "horizontal")
    grid.boxplot(unlist(ssGSEA_mat_scaled[index[i], which(group == 1)]), pos = nr-i+1 - 0.2, box_width = 0.3, 
                 gp = gpar(fill = 2), direction = "horizontal")
  }
  grid.xaxis()
  popViewport()
}
pd_status <- clinic$PFS_event[is.na(clinic$RNA_FILE)==FALSE]
names(pd_status) <- clinic$RNA_FILE[is.na(clinic$RNA_FILE)==FALSE]
hm_list <- Heatmap(
  ssGSEA_mat_scaled,
  column_split = group,
  show_column_names = TRUE,
  heatmap_legend_param = list(title = "Z-score"),
  width = unit(7, "cm"),
  top_annotation = HeatmapAnnotation(
    Group = anno_block(gp = gpar(fill = 3:2),
                       labels = c("SD/PD", "CR/PR"), 
                       labels_gp = gpar(col = "white", fontsize = 10)),
    `Disease Progression` = pd_status[colnames(ssGSEA_mat_scaled)],
    col = list(`Disease Progression` =  c(
      "0" = "#F1E290", "1" = "#F677C1"
    )),
    show_legend = FALSE
  )
) + rowAnnotation(
  boxplot = anno_multiple_boxplot,
  width = unit(3, "cm"),
  show_annotation_name = FALSE,
  show_legend = TRUE
) + rowAnnotation(
  rownames =  anno_text(rownames(ssGSEA_mat_scaled), which = "row"),
  show_annotation_name = FALSE
)

lgd_list = list(
  Legend(labels = c("SD/PD", "CR/PR"), title = "Group", 
         legend_gp = gpar(fill = 3:2)),
  Legend(labels = c("No", "Yes"), title = "Disease Progression", 
         legend_gp = gpar(fill =c("0" = "#F1E290", "1" = "#F677C1")))
)

pdf("./figures/tme_heatmap_response.pdf", height = 6, width = 10)
draw(hm_list,heatmap_legend_list = lgd_list)
dev.off()
writexl::write_xlsx(DEG_immune[[division_criteria]] %>% rownames_to_column("Cell.Type"), path = "./figures/tme_limma.xlsx")

division_criteria <- "benefit"
group <- divide_group(division_criteria)

design = model.matrix(~group)
fit = lmFit(gsva_immune, design)
fit = eBayes(fit)
DEG_immune[[division_criteria]] = topTable(fit, coef = 2, number = Inf)
head(DEG_immune[[division_criteria]],20)
tinyarray::draw_heatmap(gsva_immune[head(rownames(DEG_immune[[division_criteria]]),20),],group,show_rownames = TRUE, split_column = TRUE,show_column_title = TRUE, color = (colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlBu"))))(100), color_an = palettes(category = "box", palette = "nrc", show_col = FALSE, show_message = FALSE))

## 蛋白组学


### 蛋白组学差异分析

count_norm <- read.csv("./data/protein_matrix.csv")
conditions <- as.factor(count_norm$Group)
count_norm <- count_norm %>% column_to_rownames(var = "Name") %>% dplyr::select(-c(1:2)) %>% t() %>% as.data.frame()
#detectCores函数可以告诉你你的CPU可使用的核数
library(parallel)
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum))
clusterExport(cl = cl, varlist = c("count_norm", "conditions"))
# Run the Wilcoxon rank-sum test for each gene
pvalues <- parSapply(cl = cl, 1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=t.test(gene~conditions, data)$p.value
  rm(data)
  gc()
  return(p)
})
stopCluster(cl)
fdr=p.adjust(pvalues,method = "BH")

# Calculate fold-change for each gene

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(as.matrix(dataCon2))/rowMeans(as.matrix(dataCon1)))

# Output results based on FDR threshold

DEG<-data.frame(logFC=foldChanges, pvalue=pvalues, padj=fdr)
rownames(DEG)=rownames(count_norm)
DEG=na.omit(DEG)
DEG <- DEG[which(is.infinite(DEG$logFC) == FALSE),]
threshold<-as.factor((DEG$logFC>1|DEG$logFC<(-1)) & DEG$pvalue<0.05)

# protein_dea <- read.csv("./data/protein_dea.csv")
protein_symbol <- bitr(
  geneID = rownames(DEG),
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>% filter(UNIPROT %in% names(table(UNIPROT))[table(UNIPROT)>1] == FALSE)  %>% filter(ENTREZID %in% names(table(ENTREZID))[table(ENTREZID)>1] == FALSE)
protein_symbol_fuck <- bitr(
  geneID = rownames(DEG),
  fromType = "UNIPROT",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
) %>% filter(UNIPROT %in% names(table(UNIPROT))[table(UNIPROT)>1] == FALSE)  %>% filter(SYMBOL %in% names(table(SYMBOL))[table(SYMBOL)>1] == FALSE)
protein_dea <- DEG %>% 
  filter(rownames(.) %in% protein_symbol$UNIPROT) %>%
  # filter(protein.ID %in% names(table(protein_symbol$UNIPROT))[table(protein_symbol$UNIPROT)>1] == FALSE) %>%
  magrittr::set_rownames(protein_symbol$ENTREZID)
gene_list <- protein_dea %>%
  arrange(desc(logFC)) %>%
  dplyr::select("logFC")
gene_list_new <- gene_list[,1]
names(gene_list_new) <- rownames(gene_list)
gsea_res_pro <- list()

gsea_res_pro[["kegg2"]] <- gseKEGG(
  geneList = gene_list_new,
  organism     = "hsa",
  minGSSize    = 10,
  pvalueCutoff = 1,
  verbose      = TRUE
)

gseaplot2(gsea_res_pro[["kegg2"]], geneSetID = "hsa04612", pvalue_table = TRUE)

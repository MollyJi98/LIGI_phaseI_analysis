####################显著表型与肺癌的共定位分析####################
##step1：提取所有显著的区域（1000leadsnp），=，每个区域提取leadsnp上下游250kb的位点
##step2:在肺癌和其他肺部疾病中分别提取上述区域的位点，进行共定位分析



########################################step1########################################
R
library(data.table)
setwd('./Assoc_Lung_phenotype_low_3.8w_clean')

# 读取 extract 文件
extract <- fread("./extract_region_905_excludechrX.txt")

# 初始化变量用于存储上次读取的文件名和对应的数据
last_file <- NULL
last_data <- NULL

# 循环遍历 extract 文件的每一行,
for (i in 1:nrow(extract)) {
  # 读取当前行的信息
  file <- extract[i, FILE]
  output <- extract[i, OUTPUT]
  chr <- extract[i, CHR]
  start <- extract[i, START]
  end <- extract[i, END]
  
  # 添加输出语句，检查变量值
  cat("Processing file:", file, "\n")
  
  # 如果文件名发生变化，重新读取文件
  if (is.null(last_file) || last_file != file) {
    last_file <- file
    last_data <- fread(file)
  }
  
  # 提取符合条件的行
  data_filtered <- last_data[CHR == chr & BP >= start & BP <= end, ]
  
  # 提取指定的列
  data_filtered <- data_filtered[, c("RSID", "CHR", "BP", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", 
                                     "CHISQ_LINREG", "P_LINREG", "BETA", "SE", "CHISQ_BOLT_LMM_INF", 
                                     "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "P_BOLT_LMM"), with = FALSE]
  
  # 修改列名
  setnames(data_filtered, "RSID", "SNP")
  
  # 写入输出文件
  fwrite(data_filtered, output, sep = "\t", quote = FALSE, row.names = FALSE)
}


########################################step2########################################
############################批量写循环_OUTCOME############################
options(stringsAsFactors=F)
rm(list=ls())
library(data.table)
library(R.utils)
library(coloc)
library(dplyr)
library(stringr)
library(locuscomparer)
library(ggplot2)
library(extrafont)
library(tidyr)
library(writexl)

# 读取OUTCOME数据
OUTCOME_geno1 <- fread("outcome_dir.txt")
#读取gene数据
gene <-fread('extract_region_905_excludechrX.txt')

# 定义共定位分析函数
jc_coloc <- function(gene_id, coloc_file) {
  outcome_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_OUTCOME,
    varbeta = coloc_file$varbeta_OUTCOME,
    type = "cc",
    N = coloc_file$N,
    MAF = coloc_file$maf_OUTCOME
  )
  gwas_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_gwas,
    varbeta = coloc_file$varbeta_gwas,
    type = "quant",
    N = 35469,
    MAF = coloc_file$maf_gwas
  )
  result <- coloc.abf(dataset1 = outcome_coloc, dataset2 = gwas_coloc)
  res <- t(as.data.frame(result$summary))
  rownames(res) <- gene_id
  res <- as.data.frame(res)
  return(res)
}

# 获取所有_merge_snps.txt文件的路径
file_list <- list.files(path = "step1_output_dir",
                        pattern = "*_merge_snps.txt",
                        full.names = TRUE)

# 初始化空数据框以存储所有结果
all_res <- data.frame()

# 遍历每个文件进行分析
for (file in file_list[c(1:905)]) {
  phe <- fread(file)
  gwas_geno <- phe
  gwas_geno$chrbp = paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
  
  names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'RSID'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
  names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
  names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
  names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'
  
  gwas_geno1 <- gwas_geno %>% 
    filter((nchar(A1) == 1) & (nchar(A2) == 1)) %>%
    mutate(Base_pair = paste(A1, A2, sep = ":")) %>%
    mutate(varbeta_gwas = (SE)^2 ) %>%
    mutate(GWASsnpid = paste(CHR,BP,A2, A1, sep = ":")) %>%
    select(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas,beta_gwas, varbeta_gwas, P, RSID)
  
  gwas_geno1$Base_unified <- ifelse((gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A")),"A:C", 
                                    ifelse((gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                                           ifelse(gwas_geno1$Base_pair == "T:A","A:T",
                                                  ifelse(gwas_geno1$Base_pair == "G:C","C:G",gwas_geno1$Base_pair))))
  gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),-gwas_geno1$beta_gwas,gwas_geno1$beta_gwas)
  gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),1-gwas_geno1$maf_gwas,gwas_geno1$maf_gwas)
  gwas_geno1 <- gwas_geno1 %>%
    mutate(SNP = paste(chrbp, Base_unified, sep = ":"))
  
  # 合并GWAS数据
  coloc_data <- OUTCOME_geno1 %>% inner_join(gwas_geno1, by = "SNP")

  #匹配gene
  basename<- basename(file)
  basename<-sub("_merge_snps.txt$", "", basename)
  # 分割 basename_no_ext
  parts <- unlist(strsplit(basename, "\\."))
  # 提取 phe, CHR, BP
  phe <- paste0(parts[1], "_assoc.result")
  CHR_value <- as.numeric(parts[2])
  BP_value <- as.numeric(parts[3])

  # 在 gene 数据框中查找匹配的行
  matching_row <- subset(gene,gene$FILE == phe & gene$CHR == CHR_value & gene$BP == BP_value)
  # 提取 GENE 列内容
  GENE <- ifelse(nrow(matching_row) > 0, matching_row$GENE, NA)

  # 共定位分析
  res <- jc_coloc(GENE, coloc_data) 
  
  # 合并结果
  all_res <- rbind(all_res, res)
  
  # 绘制区域图
  outcome_coloc <- coloc_data
  
  outcome_fn <- outcome_coloc[, c('RSID', 'P_OUTCOME')]
  names(outcome_fn) <- c('rsid', 'pval')
  gwas_fn <- outcome_coloc[, c('RSID', 'P')]
  names(gwas_fn) <- c('rsid', 'pval')
  
  # 找最显著的点
  data <- cbind(outcome_fn, gwas_fn[, 2])
  names(data) <- c('RSID', 'P_OUTCOME', 'P_phe')
  data$sum_p <- (-log(data$P_OUTCOME)) + (-log(data$P_phe))
  rsid_max_sum_p <- data[which.max(data$sum_p), "RSID"]
  
  setwd('outputdir') #set output dir
  locus_plot = paste0("coloc_",GENE,"_",basename, ".pdf")
  marker_col = "rsid"
  pval_col = "pval"
  leadsnp = rsid_max_sum_p
  
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = outcome_fn, 
                    population = "EAS", 
                    title1 = "Phenotype", title2 = "OUTCOME", 
                    legend_position = 'topleft', 
                    marker_col1 = marker_col, pval_col1 = pval_col, marker_col2 = marker_col, pval_col2 = pval_col,
                    snp = leadsnp,
                    genome = "hg19")
  p <- p + theme(
    text = element_text(family = "Arial"),  # 设置整体字体为 Arial
    plot.title = element_text(size = 15, family = "Arial"),  # 设置标题字体为 Arial
    axis.title = element_text(size = 12, family = "Arial"),  # 设置轴标题字体为 Arial
    axis.text = element_text(size = 10, family = "Arial"),   # 设置轴文本字体为 Arial
    legend.text = element_text(size = 10, family = "Arial"),  # 设置图例文本字体为 Arial
    legend.title = element_text(size = 12, family = "Arial")  # 设置图例标题字体为 Arial
  )
  ggsave(locus_plot, plot = p, width = 9, height = 5, device = cairo_pdf)
}

all_res$pheno<-file_list
# 将合并的结果保存到文件
write.table(all_res, "all_coloc_results_OUTCOME.txt", quote = F, sep = "\t", row.names = T)





###########整理结果###########
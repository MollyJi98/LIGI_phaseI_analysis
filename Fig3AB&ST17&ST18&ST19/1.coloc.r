cd ./Assoc_Lung_phenotype_coloc

#####R
options(stringsAsFactors=F)
rm(list=ls())
library(data.table)
library(tidyfst)
library(R.utils)
library(coloc)


filenames <- list.files(path="./Assoc_Lung_phenotype_low_3.8w_clean", pattern="*_assoc.result") 
file=data.frame(filenames)
file$names <- gsub( "_assoc.result", "",file$filenames)
file=data.table(file)
dim(file) #160
file$N=338
write.table(file,"./filelist_clean.list",quote = F,sep = "\t",row.names = F,col.names = T)



#############################################################################data prepare demo
library(data.table)
library(tidyfst)
library(R.utils)
library(coloc)

eqtl_data <- parse_fst("data/338sample_bulkeqtl_egene_data.fst") 
list <- fread("filelist_clean.list")
gwasdir = 'data/Assoc_Lung_phenotype_low_3.8w_clean'  
datadir = 'data/Assoc_Lung_phenotype_coloc/data'

# 举例：选择一个表型进行分析
phenoname <- list$names[1] 

# 读取
gwas_geno <- fread(sprintf("%s/%s_assoc.result", gwasdir, phenoname)) 
head(gwas_geno)

# 处理变量名
gwas_geno$chrbp = paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'GWASsnpid'
names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'

gwas_geno1 <- gwas_geno %>% 
  filter_dt((nchar(A1) == 1) & (nchar(A2) == 1)) %>%
  mutate_dt(Base_pair = paste(A1, A2, sep = ":")) %>%
  mutate_dt(varbeta_gwas = (SE)^2 ) %>%
  select_mix(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas, beta_gwas, varbeta_gwas, P)

# 处理等位基因
gwas_geno1$Base_unified <- ifelse((gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A")), "A:C", 
                              ifelse((gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                              ifelse(gwas_geno1$Base_pair == "T:A", "A:T",
                              ifelse(gwas_geno1$Base_pair == "G:C", "C:G", gwas_geno1$Base_pair))))

gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 
                                -gwas_geno1$beta_gwas, gwas_geno1$beta_gwas)
gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 
                               1-gwas_geno1$maf_gwas, gwas_geno1$maf_gwas)

gwas_geno1 <- gwas_geno1 %>% mutate_dt(SNP = paste(chrbp, Base_unified, sep = ":"))

# 合并GWAS和eQTL数据
coloc_data <- eqtl_data %>% inner_join_dt(gwas_geno1, by = "SNP")  # 保留GWAS-eQTL共有的SNP

# 导出结果
export_fst(coloc_data, sprintf("%s/%s.fst", datadir, phenoname))

# 清理内存
rm(gwas_geno, gwas_geno1, coloc_data)




##############################################coloc 示例代码##############################################
library(data.table)
library(tidyfst)
library(coloc)
library(tibble)
library(dplyr)

# 设置路径
datadir <- 'data/Assoc_Lung_phenotype_coloc/data'
resultdir <- 'data/Assoc_Lung_phenotype_coloc/result'
plotdir <- 'data/Assoc_Lung_phenotype_coloc/plot'

# 定义共定位分析函数
qtl_coloc <- function(gene_id, coloc_file) {
  eqtl_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_eqtl,
    varbeta = coloc_file$varbeta_eqtl,
    type = "quant",
    N = 338,
    MAF = coloc_file$maf
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
  result <- coloc.abf(dataset1 = eqtl_coloc, dataset2 = gwas_coloc)
  res <- t(as.data.frame(result$summary))
  rownames(res) <- gene_id
  res <- as.data.frame(res)
  return(res)
}

# 读取并整理数据
list <- fread('data/filelist_clean.list')
phenoname <- list$names[1]  # 选择第一个表型进行分析

dat <- parse_fst(sprintf("%s/%s.fst", datadir, phenoname)) %>%
  as.data.frame() %>%
  arrange_dt(gene_id)

# 基因列表
gene_list <- unique(dat$gene_id)
res_all <- data.frame()

# 对每个基因执行共定位分析
for (gene_id in gene_list) {
  sub <- dat %>%
    filter_dt(gene_id == gene_id) %>%
    distinct()  # 删除重复行
  if (nrow(sub) > 30) {
    res <- qtl_coloc(gene_id, sub)
    res_all <- rbind(res_all, res)
  }
}

# 保存共定位分析结果
res_all <- res_all %>%
  rownames_to_column(var = "gene_id")
fwrite(res_all, sprintf("%s/%s_colocresult.csv", resultdir, phenoname))

# 对所有表型进行共定位分析后的可视化
egenecoloc2 <- fread(sprintf("%s/plot/allegene_coloc_sigresult.csv", resultdir))
egenecoloc_sig <- subset(egenecoloc2, PPH4 >= 0.7)
fwrite(egenecoloc_sig, sprintf("%s/allegene_coloc_sigresult.csv", plotdir), sep = ",", quote = FALSE)

# 绘制区域图
library(locuscomparer)
for (gene_name in unique(egenecoloc_sig$gene_name)) {
  myphenotype <- list$names[1]
  subdat <- siggene_allsnp %>%
    filter_dt(gene_name == gene_name)
  
  gwas_fn <- sprintf("%s/%s/%s_GWAS.txt", plotdir, myphenotype, gene_name)
  eqtl_fn <- sprintf("%s/%s/%s_eQTL.txt", plotdir, myphenotype, gene_name)
  
  locus_plot <- sprintf("%s/%s/%s_coloc.png", plotdir, myphenotype, gene_name)
  
  fwrite(gwas, gwas_fn, sep = "\t", quote = FALSE)
  fwrite(eqtl, eqtl_fn, sep = "\t", quote = FALSE)
  
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, population = "EAS", title1 = "GWAS", title2 = "eQTL",
                    legend_position = 'topleft', marker_col1 = "rsid", pval_col1 = "pval", marker_col2 = "rsid",
                    pval_col2 = "pval", snp = gene_name, genome = "hg19")
  
  ggsave(locus_plot, p, width = 9, height = 5)
}

# 最终汇总
df2 <- data.table(df2)
df2[, "EA" := tstrsplit(SNP, ":")[3]]
df2[, "OA" := tstrsplit(SNP, ":")[4]]
df2$OR <- round(exp(df2$beta_gwas), 2)
df2$SEbeta <- sqrt(df2$varbeta_gwas)
df2$L95 <- round(exp(df2$beta_gwas - df2$SEbeta * 1.96), 2)
df2$U95 <- round(exp(df2$beta_gwas + df2$SEbeta * 1.96), 2)
df2$OR_GWAS <- paste0(df2$OR, "(", df2$L95, ",", df2$U95, ")")
df3 <- df2[, c("phenotype", "gene_name", "gene_id", "seqnames", "Cytoband", "nsnps", "PPH3", "PPH4", "rsid", "CHR", "BP", "EA", "OA", "maf_gwas", "OR_GWAS", "P", "beta_eqtl", "pval_nominal")]

# 保存最终汇总
fwrite(df3, sprintf("%s/plot/allLungphenotype_colocsig_topsnp_log1.xls", plotdir), sep = "\t", quote = FALSE)

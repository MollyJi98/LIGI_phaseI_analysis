##########已经整理成位点-基因-最显著的影像组特征的表单，写循环进行分析#############

library(data.table)
library(R.utils)
library(coloc)
library(dplyr)
library(stringr)
library(locuscomparer)
library(ggplot2)
library(extrafont)

# define function
jc_coloc <- function(gene_id, coloc_file) {
  outcome_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_NSCLC,
    varbeta = coloc_file$varbeta_NSCLC,
    type = "cc",
    N = coloc_file$N,
    MAF = coloc_file$maf_NSCLC
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

result<-data.frame() #储存共定位结果
result2<-data.frame() #储存最显著的snp
NSCLC_geno1<-fread('3_China_pop_NSCLC_updated_0916_for_coloc.txt')


# 遍历 min_p_bolt_lmm_data 的每一行
for (i in 1:nrow(min_p_bolt_lmm_data)) {
  
  # 获取当前行的信息
  current_row <- min_p_bolt_lmm_data[i, ]
  filename <- current_row$filename
  chr_value <- current_row$CHR
  bp_value <- current_row$BP
  gene_value <- current_row$gene 
  
  # 生成文件路径
  file_path <- file.path("./Assoc_Lung_phenotype_low_3.8w_clean", filename)
  
  # 读取文件
  gwas_geno <- fread(file_path)
  gwas_geno$chrbp <- paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
  
  # 重命名列
  names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'GWASsnpid'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
  names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
  names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
  names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'
  
  # 过滤和转换数据
  gwas_geno1 <- gwas_geno %>%
    filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
    mutate(Base_pair = paste(A1, A2, sep = ":")) %>%
    mutate(varbeta_gwas = (SE)^2) %>%
    select(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas, beta_gwas, varbeta_gwas, P, RSID)
  
  # 整理等位基因
  gwas_geno1$Base_unified <- ifelse(gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A"), "A:C", 
                                    ifelse(gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A"), "A:G",
                                    ifelse(gwas_geno1$Base_pair == "T:A", "A:T",
                                    ifelse(gwas_geno1$Base_pair == "G:C", "C:G", gwas_geno1$Base_pair))))
  gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), -gwas_geno1$beta_gwas, gwas_geno1$beta_gwas)
  gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 1 - gwas_geno1$maf_gwas, gwas_geno1$maf_gwas)
  gwas_geno1 <- gwas_geno1 %>%
    mutate(SNP = paste(chrbp, Base_unified, sep = ":"))
  
  # 合并GWAS数据
  coloc_data <- NSCLC_geno1 %>% inner_join(gwas_geno1, by = "SNP")
  
  # 筛选目标区域数据
  coloc_file2 <- coloc_data %>%
    filter(CHR == chr_value & BP > bp_value - 500000 & BP < bp_value + 500000)
  
  # 进行共定位分析
  res <- jc_coloc(gene_value, coloc_file2)
  result <- rbind(result, res)
  
  # 画区域图
  outcome_coloc <- coloc_file2
  outcome_fn <- outcome_coloc %>% select(RSID, P_NSCLC)
  names(outcome_fn) <- c('rsid', 'pval')
  gwas_fn <- outcome_coloc %>% select(RSID, P)
  names(gwas_fn) <- c('rsid', 'pval')
  
  # 找最显著的点
  data <- cbind(outcome_fn, gwas_fn[, 2])
  names(data) <- c('RSID', 'P_NSCLC', 'P_phe')
  data$sum_p <- -log(data$P_NSCLC) + -log(data$P_phe)
  rsid_max_sum_p <- data[which.max(data$sum_p), "RSID"]
  result2<-rbind(result2,rsid_max_sum_p)
  
  locus_plot <- paste0(gene_value, "_lungcancer_", gsub("_assoc.result", "", filename), "_coloc.pdf")
  
  # 画图
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = outcome_fn, 
                    population = "EAS", 
                    title1 = "Phenotype", title2 = "NSCLC", 
                    legend_position = 'topleft', 
                    marker_col1 = marker_col, pval_col1 = pval_col, marker_col2 = marker_col, pval_col2 = pval_col,
                    snp = rsid_max_sum_p,
                    genome = "hg19")+ theme(
    text = element_text(family = "Arial"),  
    plot.title = element_text(size = 15, family = "Arial"),  
    axis.title = element_text(size = 12, family = "Arial"), 
    axis.text = element_text(size = 10, family = "Arial"),  
    legend.text = element_text(size = 10, family = "Arial"), 
    legend.title = element_text(size = 12, family = "Arial") 
  )
  
  # 保存图像
  ggsave(locus_plot, plot = p, width = 9, height = 5, device = cairo_pdf)
  
}

# 输出最终结果
print(result)
print(result2)
result<-cbind(result,result2)
result$gene<-rownames(result)
result<-result[,c(8,1:7)]

write_xlsx(result, "result_14loci_500kb_coloc.xlsx")








####################################################COPD####################################################
R
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

############################整理COPD数据############################
COPD <-fread("/data1/Imaging_Assoc_2023/calculation_result_JC/LDSC/lung_related_trait/copd_meta/copd_meta_forLDSC.txt") 
dim(COPD) #10856597         10

COPD_geno<-COPD
head(COPD_geno)
COPD_geno$chrbp = paste(COPD_geno$CHR, COPD_geno$BP, sep = ":")

# gwas_geno=gwas_geno[!gwas_geno$chrbp %in% removesnp$CHR_BP ,] #24,470,237 SNPs
names(COPD_geno)[names(COPD_geno) == 'SNP'] <- 'COPDsnpid'
names(COPD_geno)[names(COPD_geno) == 'FRQ'] <- 'maf_COPD'
names(COPD_geno)[names(COPD_geno) == 'BETA'] <- 'beta_COPD'
names(COPD_geno)[names(COPD_geno) == 'P'] <- 'P_COPD'

# Base_pair 为A1:A2 (effect allele:ref allele)
COPD_geno1 <- COPD_geno %>%
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  mutate(Base_pair = paste(A1, A2, sep = ":"),
         varbeta_COPD = SE^2) %>%
  select(CHR, BP, COPDsnpid, chrbp, Base_pair, maf_COPD, beta_COPD, varbeta_COPD, P_COPD,N)

# dim(COPD_geno1)
rm(COPD_geno)

# 整理等位基因
COPD_geno1$Base_unified <- ifelse((COPD_geno1$Base_pair %in% c("T:G", "G:T", "C:A")),"A:C", 
                             ifelse((COPD_geno1$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                             ifelse(COPD_geno1$Base_pair == "T:A","A:T",
                             ifelse(COPD_geno1$Base_pair == "G:C","C:G",COPD_geno1$Base_pair))))
COPD_geno1$beta_COPD <- ifelse(COPD_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),-COPD_geno1$beta_COPD,COPD_geno1$beta_COPD)
COPD_geno1$maf_COPD <- ifelse(COPD_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),1-COPD_geno1$maf_COPD,COPD_geno1$maf_COPD)
COPD_geno1 <- COPD_geno1 %>%
  mutate(SNP = paste(chrbp, Base_unified, sep = ":")) %>%
  select(SNP, COPDsnpid, Base_pair, maf_COPD, beta_COPD, varbeta_COPD, P_COPD,N)

names(COPD_geno1)[3]<-"Base_pair_COPD"

write.table(COPD_geno1,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_for_coloc.txt",quote=F,sep="\t",row.names=F)


############################批量写循环############################
#提取COPD显著的位点
data<-subset(COPD,COPD$P<5e-8)
rsid<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/5.LDSC/step2_outcome_prepare/COPD_for_LDSC.txt')
rsid$base = paste(toupper(rsid$A1),toupper(rsid$A2),sep=':')
rsid$base_1<-ifelse((rsid$base=='T:G'|rsid$base=='G:T'|rsid$base=='C:A'),'A:C',ifelse((rsid$base=='T:C'|rsid$base=='C:T'|rsid$base=='G:A'),'A:G',ifelse(rsid$base=='T:A','A:T',ifelse(rsid$base=='G:C','C:G',rsid$base))))
rsid$CHR_BP_ALLELE =paste(rsid$CHR,rsid$BP,rsid$base_1,sep=':')
rsid<-subset(rsid,rsid$CHR_BP_ALLELE %in% data$SNP)

data1<-merge(data,rsid[,c("CHR_BP_ALLELE","SNP")],by.x="SNP",by.y="CHR_BP_ALLELE")
names(data1)[names(data1) == 'SNP'] <- 'suibian'
names(data1)[names(data1) == 'SNP.y'] <- 'SNP'

write.table(data1[,c("SNP","P")],'/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_sig_snps.txt',sep="\t",quote=F,row.names=F)

#clumping
plink --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS \
      --clump /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_sig_snps.txt \
      --clump-kb 1000 \
      --clump-p1 0.05 \
      --clump-r2 0.1 \
      --out /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_sig_snps
#一共10个lead snps#



############################整理ILD数据############################
ILD <-fread("/data1/Imaging_Assoc_2023/calculation_result_JC/LDSC/lung_related_trait/ILD_meta/ILD_meta_forLDSC.txt") 
dim(ILD) #10856597         10

ILD_geno<-ILD
head(ILD_geno)
ILD_geno$chrbp = paste(ILD_geno$CHR, ILD_geno$BP, sep = ":")

# gwas_geno=gwas_geno[!gwas_geno$chrbp %in% removesnp$CHR_BP ,] #24,470,237 SNPs
names(ILD_geno)[names(ILD_geno) == 'SNP'] <- 'ILDsnpid'
names(ILD_geno)[names(ILD_geno) == 'FRQ'] <- 'maf_ILD'
names(ILD_geno)[names(ILD_geno) == 'BETA'] <- 'beta_ILD'
names(ILD_geno)[names(ILD_geno) == 'P'] <- 'P_ILD'

# Base_pair 为A1:A2 (effect allele:ref allele)
ILD_geno1 <- ILD_geno %>%
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  mutate(Base_pair = paste(A1, A2, sep = ":"),
         varbeta_ILD = SE^2) %>%
  select(CHR, BP, ILDsnpid, chrbp, Base_pair, maf_ILD, beta_ILD, varbeta_ILD, P_ILD,N)

# dim(ILD_geno1)
rm(ILD_geno)

# 整理等位基因
ILD_geno1$Base_unified <- ifelse((ILD_geno1$Base_pair %in% c("T:G", "G:T", "C:A")),"A:C", 
                             ifelse((ILD_geno1$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                             ifelse(ILD_geno1$Base_pair == "T:A","A:T",
                             ifelse(ILD_geno1$Base_pair == "G:C","C:G",ILD_geno1$Base_pair))))
ILD_geno1$beta_ILD <- ifelse(ILD_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),-ILD_geno1$beta_ILD,ILD_geno1$beta_ILD)
ILD_geno1$maf_ILD <- ifelse(ILD_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),1-ILD_geno1$maf_ILD,ILD_geno1$maf_ILD)
ILD_geno1 <- ILD_geno1 %>%
  mutate(SNP = paste(chrbp, Base_unified, sep = ":")) %>%
  select(SNP, ILDsnpid, Base_pair, maf_ILD, beta_ILD, varbeta_ILD, P_ILD,N)

names(ILD_geno1)[3]<-"Base_pair_ILD"

write.table(ILD_geno1,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_for_coloc.txt",quote=F,sep="\t",row.names=F)

############################批量写循环############################
#提取ILD显著的位点
data<-subset(ILD,ILD$P<5e-8)
rsid<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/5.LDSC/step2_outcome_prepare/ILD_for_LDSC.txt')
rsid$base = paste(toupper(rsid$A1),toupper(rsid$A2),sep=':')
rsid$base_1<-ifelse((rsid$base=='T:G'|rsid$base=='G:T'|rsid$base=='C:A'),'A:C',ifelse((rsid$base=='T:C'|rsid$base=='C:T'|rsid$base=='G:A'),'A:G',ifelse(rsid$base=='T:A','A:T',ifelse(rsid$base=='G:C','C:G',rsid$base))))
rsid$CHR_BP_ALLELE =paste(rsid$CHR,rsid$BP,rsid$base_1,sep=':')
rsid<-subset(rsid,rsid$CHR_BP_ALLELE %in% data$SNP)

data1<-merge(data,rsid[,c("CHR_BP_ALLELE","SNP")],by.x="SNP",by.y="CHR_BP_ALLELE")
names(data1)[names(data1) == 'SNP'] <- 'suibian'
names(data1)[names(data1) == 'SNP.y'] <- 'SNP'

write.table(data1[,c("SNP","P")],'/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_sig_snps.txt',sep="\t",quote=F,row.names=F)

#clumping
plink --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS \
      --clump /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_sig_snps.txt \
      --clump-kb 1000 \
      --clump-p1 0.05 \
      --clump-r2 0.1 \
      --out /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_sig_snps
#一共2个lead snps#





############################整理ASTHMA数据############################
ASTHMA <-fread("/data1/Imaging_Assoc_2023/calculation_result_JC/LDSC/lung_related_trait/Asthma_meta/Asthma_meta_forLDSC.txt") 
dim(ASTHMA) #10856597         10

ASTHMA_geno<-ASTHMA
head(ASTHMA_geno)
ASTHMA_geno$chrbp = paste(ASTHMA_geno$CHR, ASTHMA_geno$BP, sep = ":")

# gwas_geno=gwas_geno[!gwas_geno$chrbp %in% removesnp$CHR_BP ,] #24,470,237 SNPs
names(ASTHMA_geno)[names(ASTHMA_geno) == 'SNP'] <- 'ASTHMAsnpid'
names(ASTHMA_geno)[names(ASTHMA_geno) == 'FRQ'] <- 'maf_ASTHMA'
names(ASTHMA_geno)[names(ASTHMA_geno) == 'BETA'] <- 'beta_ASTHMA'
names(ASTHMA_geno)[names(ASTHMA_geno) == 'P'] <- 'P_ASTHMA'

# Base_pair 为A1:A2 (effect allele:ref allele)
ASTHMA_geno1 <- ASTHMA_geno %>%
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  mutate(Base_pair = paste(A1, A2, sep = ":"),
         varbeta_ASTHMA = SE^2) %>%
  select(CHR, BP, ASTHMAsnpid, chrbp, Base_pair, maf_ASTHMA, beta_ASTHMA, varbeta_ASTHMA, P_ASTHMA,N)

# dim(ASTHMA_geno1)
rm(ASTHMA_geno)

# 整理等位基因
ASTHMA_geno1$Base_unified <- ifelse((ASTHMA_geno1$Base_pair %in% c("T:G", "G:T", "C:A")),"A:C", 
                             ifelse((ASTHMA_geno1$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                             ifelse(ASTHMA_geno1$Base_pair == "T:A","A:T",
                             ifelse(ASTHMA_geno1$Base_pair == "G:C","C:G",ASTHMA_geno1$Base_pair))))
ASTHMA_geno1$beta_ASTHMA <- ifelse(ASTHMA_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),-ASTHMA_geno1$beta_ASTHMA,ASTHMA_geno1$beta_ASTHMA)
ASTHMA_geno1$maf_ASTHMA <- ifelse(ASTHMA_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"),1-ASTHMA_geno1$maf_ASTHMA,ASTHMA_geno1$maf_ASTHMA)
ASTHMA_geno1 <- ASTHMA_geno1 %>%
  mutate(SNP = paste(chrbp, Base_unified, sep = ":")) %>%
  select(SNP, ASTHMAsnpid, Base_pair, maf_ASTHMA, beta_ASTHMA, varbeta_ASTHMA, P_ASTHMA,N)

names(ASTHMA_geno1)[3]<-"Base_pair_ASTHMA"

write.table(ASTHMA_geno1,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_for_coloc.txt",quote=F,sep="\t",row.names=F)

############################批量写循环############################
#提取ASTHMA显著的位点
data<-subset(ASTHMA,ASTHMA$P<5e-8)
rsid<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/5.LDSC/step2_outcome_prepare/Asthma_for_LDSC.txt')
rsid$base = paste(toupper(rsid$A1),toupper(rsid$A2),sep=':')
rsid$base_1<-ifelse((rsid$base=='T:G'|rsid$base=='G:T'|rsid$base=='C:A'),'A:C',ifelse((rsid$base=='T:C'|rsid$base=='C:T'|rsid$base=='G:A'),'A:G',ifelse(rsid$base=='T:A','A:T',ifelse(rsid$base=='G:C','C:G',rsid$base))))
rsid$CHR_BP_ALLELE =paste(rsid$CHR,rsid$BP,rsid$base_1,sep=':')
rsid<-subset(rsid,rsid$CHR_BP_ALLELE %in% data$SNP)

data1<-merge(data,rsid[,c("CHR_BP_ALLELE","SNP")],by.x="SNP",by.y="CHR_BP_ALLELE")
names(data1)[names(data1) == 'SNP'] <- 'suibian'
names(data1)[names(data1) == 'SNP.y'] <- 'SNP'

write.table(data1[,c("SNP","P")],'/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_sig_snps.txt',sep="\t",quote=F,row.names=F)

#clumping
plink --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS \
      --clump /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_sig_snps.txt \
      --clump-kb 1000 \
      --clump-p1 0.05 \
      --clump-r2 0.1 \
      --out /data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_sig_snps
#一共10个lead snps#





#对于clumping后的lead snp,在160个表型中提取10+2+76个点的p值，找到在哪个表型里最显著

# 定义源文件夹和目标文件夹  
source_folder <- "/data1/Imaging_Assoc_2023/Assoc_Lung_phenotype_low_3.8w_clean"  
target_folder <- "/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/data/COPD_ILD_ASTH"  
  
# 读取要提的点
data1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
data2<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
data3<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
data1<-rbind(data1,data2)
data1<-rbind(data1,data3)
data1$chr_bp<-paste(data1$CHR,data1$BP,sep=":") 
write.table(data1,'/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/data/COPD_ILD_ASTH/3_outcome_88_sig_snps.txt',sep="\t",quote=F,row.names=F)


# 遍历源文件夹下的所有文件  
files <- list.files(path = source_folder, pattern = "\\.result$", full.names = TRUE)  
  
for (file in files) {  
  # 读取当前文件  
  assoc_df <- fread(file)  
    
  # 提取文件名（不带_assoc.result）  
  base_filename <- gsub("_assoc.result", "", basename(file))  
    
  # 合并文件名和输出路径  
  output_file <- file.path(target_folder, paste0(base_filename, "_88_leadSNPS.txt"))  
    
  # 提取与GWAS catalog匹配的chrpos  
  matched_rows <- subset(assoc_df,assoc_df$chrpos %in% data1$chr_bp)
    
  # 提取chrpos, BETA, SE, P_BOLT_LMM列，并保存到文件  
  write.table(matched_rows[, c("chrpos", "BETA", "SE", "P_BOLT_LMM","RSID")],   
              file = output_file,   
              sep = "\t",   
              quote = FALSE,   
              row.names = FALSE)  
}  
  
# 脚本结束  
print("All files have been processed.")

#读取并合并到一个大数据框中#
file_path <- "/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/data/COPD_ILD_ASTH"

# 获取所有文件名
file_names <- list.files(file_path, pattern = "_88_leadSNPS.txt$", full.names = TRUE)

# 初始化一个空列表，用于存储每个文件的数据框
data_list <- list()

# 读取每个文件并处理
for (file in file_names) {
  # 读取文件
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 提取文件名并修改
  filename <- basename(file)
  filename <- sub("_88_leadSNPS.txt$", "_assoc.result", filename)
  
  # 添加新的列
  df$filename <- filename
  
  # 将数据框添加到列表中
  data_list[[length(data_list) + 1]] <- df
}

# 将所有数据框合并成一个大数据框
combined_data <- do.call(rbind, data_list)

# 转换 P_BOLT_LMM 列为数值类型，以确保可以进行比较
combined_data$P_BOLT_LMM <- as.numeric(combined_data$P_BOLT_LMM)

# 按 RSID 分组，找到每个分组中 P_BOLT_LMM 最小的那一行
min_p_bolt_lmm_data <- combined_data %>%
  group_by(RSID) %>%
  filter(P_BOLT_LMM == min(P_BOLT_LMM)) %>%
  select(RSID, P_BOLT_LMM, filename, chrpos) %>%
  separate(chrpos, into = c("CHR", "BP"), sep = ":", convert = TRUE) %>%
  select(RSID, P_BOLT_LMM, filename, CHR, BP)  %>% # 去掉 chrpos 列
  as.data.frame()


#查找gene（使用vep在线版查找了88个位点的对应的gene）
gene<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/data/3_outcome_88_sig_snps_gene.txt')[,c("SNP","gene")]
min_p_bolt_lmm_data<-merge(min_p_bolt_lmm_data,gene,by.x="RSID",by.y="SNP")
head(min_p_bolt_lmm_data) 
table(duplicated(min_p_bolt_lmm_data$RSID)) #一共81个点在我们的数据库中,且min_p_bolt_lmm_data有的行重复（p一样）


##########已经整理成位点-基因-最显著的影像组特征的表单，接下来要批量写循环！！！！#############
library(data.table)
library(R.utils)
library(coloc)
library(dplyr)
library(stringr)
library(locuscomparer)
library(ggplot2)
library(extrafont)

# define function
jc_coloc <- function(gene_id, coloc_file) {
  outcome_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_COPD,
    varbeta = coloc_file$varbeta_COPD,
    type = "cc",
    N = coloc_file$N,
    MAF = coloc_file$maf_COPD
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

result<-data.frame() #储存共定位结果
result2<-data.frame() #储存最显著的snp

############################################COPD############################################
data1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
COPD_geno1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/COPD_META_for_coloc.txt')
min_p_bolt_lmm_data_COPD<-subset(min_p_bolt_lmm_data,min_p_bolt_lmm_data$RSID %in% data1$SNP)

write.table(min_p_bolt_lmm_data_COPD,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/data/COPD_MATCH_phenotype.txt",quote=F,sep="\t",row.names=F)

# 设置图形保存路径
setwd('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/data')
# 遍历 min_p_bolt_lmm_data_COPD 的1:10行
for (i in 3:nrow(min_p_bolt_lmm_data_COPD)) {
  
  # 获取当前行的信息
  current_row <- min_p_bolt_lmm_data_COPD[i, ]
  filename <- current_row$filename
  chr_value <- current_row$CHR
  bp_value <- current_row$BP
  gene_value <- current_row$gene 
  
  # 生成文件路径
  file_path <- file.path("/data1/Imaging_Assoc_2023/Assoc_Lung_phenotype_low_3.8w_clean", filename)
  
  # 读取文件
  gwas_geno <- fread(file_path)
  gwas_geno$chrbp <- paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
  
  # 重命名列
  names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'GWASsnpid'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
  names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
  names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
  names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'
  
  # 过滤和转换数据
  gwas_geno1 <- gwas_geno %>%
    filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
    mutate(Base_pair = paste(A1, A2, sep = ":")) %>%
    mutate(varbeta_gwas = (SE)^2) %>%
    select(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas, beta_gwas, varbeta_gwas, P, RSID)
  
  # 整理等位基因
  gwas_geno1$Base_unified <- ifelse(gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A"), "A:C", 
                                    ifelse(gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A"), "A:G",
                                    ifelse(gwas_geno1$Base_pair == "T:A", "A:T",
                                    ifelse(gwas_geno1$Base_pair == "G:C", "C:G", gwas_geno1$Base_pair))))
  gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), -gwas_geno1$beta_gwas, gwas_geno1$beta_gwas)
  gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 1 - gwas_geno1$maf_gwas, gwas_geno1$maf_gwas)
  gwas_geno1 <- gwas_geno1 %>%
    mutate(SNP = paste(chrbp, Base_unified, sep = ":"))
  
  # 合并GWAS数据
  coloc_data <- COPD_geno1 %>% inner_join(gwas_geno1, by = "SNP")
  
  # 筛选目标区域数据
  coloc_file2 <- coloc_data %>%
    filter(CHR == chr_value & BP > bp_value - 500000 & BP < bp_value + 500000)
  
  # 进行共定位分析
  res <- jc_coloc(gene_value, coloc_file2)
  result <- rbind(result, res)
  
  # 画区域图
  outcome_coloc <- coloc_file2
  outcome_fn <- outcome_coloc %>% select(RSID, P_COPD)
  names(outcome_fn) <- c('rsid', 'pval')
  gwas_fn <- outcome_coloc %>% select(RSID, P)
  names(gwas_fn) <- c('rsid', 'pval')
  
  # 找最显著的点
  data <- cbind(outcome_fn, gwas_fn[, 2])
  names(data) <- c('RSID', 'P_COPD', 'P_phe')
  data$sum_p <- -log(data$P_COPD) + -log(data$P_phe)
  rsid_max_sum_p <- data[which.max(data$sum_p), "RSID"]
  result2<-rbind(result2,rsid_max_sum_p)
  
  locus_plot <- paste0(gene_value, "_COPD_", gsub("_assoc.result", "", filename), "_coloc.pdf")
  
  # 画图
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = outcome_fn, 
                    population = "EAS", 
                    title1 = "Phenotype", title2 = "COPD", 
                    legend_position = 'topleft', 
                    marker_col1 = marker_col, pval_col1 = pval_col, marker_col2 = marker_col, pval_col2 = pval_col,
                    snp = rsid_max_sum_p,
                    genome = "hg19")+ theme(
    text = element_text(family = "Arial"),  
    plot.title = element_text(size = 15, family = "Arial"),  
    axis.title = element_text(size = 12, family = "Arial"), 
    axis.text = element_text(size = 10, family = "Arial"),  
    legend.text = element_text(size = 10, family = "Arial"), 
    legend.title = element_text(size = 12, family = "Arial") 
  )
  
  # 保存图像
  ggsave(locus_plot, plot = p, width = 9, height = 5, device = cairo_pdf)
  
}

# 输出最终结果
print(result)
print(result2)
result<-cbind(result,result2)
result$gene<-rownames(result)
result<-result[,c(8,1:7)]

write_xlsx(result, "/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/COPD/data/result_10loci_500kb_coloc.xlsx")






############################################ASTHMA############################################
# define function
jc_coloc <- function(gene_id, coloc_file) {
  outcome_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_ASTHMA,
    varbeta = coloc_file$varbeta_ASTHMA,
    type = "cc",
    N = coloc_file$N,
    MAF = coloc_file$maf_ASTHMA
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

result<-data.frame() #储存共定位结果
result2<-data.frame() #储存最显著的snp


data1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
ASTHMA_geno1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/ASTHMA_META_for_coloc.txt')
min_p_bolt_lmm_data_ASTHMA<-subset(min_p_bolt_lmm_data,min_p_bolt_lmm_data$RSID %in% data1$SNP)

write.table(min_p_bolt_lmm_data_ASTHMA,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/data/ASTHMA_MATCH_phenotype.txt",quote=F,sep="\t",row.names=F)

# 设置图形保存路径
setwd('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/data')
# 遍历 min_p_bolt_lmm_data_ASTHMA 
for (i in 1:nrow(min_p_bolt_lmm_data_ASTHMA)) {
  
  # 获取当前行的信息
  current_row <- min_p_bolt_lmm_data_ASTHMA[i, ]
  filename <- current_row$filename
  chr_value <- current_row$CHR
  bp_value <- current_row$BP
  gene_value <- current_row$gene 
  
  # 生成文件路径
  file_path <- file.path("/data1/Imaging_Assoc_2023/Assoc_Lung_phenotype_low_3.8w_clean", filename)
  
  # 读取文件
  gwas_geno <- fread(file_path)
  gwas_geno$chrbp <- paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
  
  # 重命名列
  names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'GWASsnpid'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
  names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
  names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
  names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'
  
  # 过滤和转换数据
  gwas_geno1 <- gwas_geno %>%
    filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
    mutate(Base_pair = paste(A1, A2, sep = ":")) %>%
    mutate(varbeta_gwas = (SE)^2) %>%
    select(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas, beta_gwas, varbeta_gwas, P, RSID)
  
  # 整理等位基因
  gwas_geno1$Base_unified <- ifelse(gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A"), "A:C", 
                                    ifelse(gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A"), "A:G",
                                    ifelse(gwas_geno1$Base_pair == "T:A", "A:T",
                                    ifelse(gwas_geno1$Base_pair == "G:C", "C:G", gwas_geno1$Base_pair))))
  gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), -gwas_geno1$beta_gwas, gwas_geno1$beta_gwas)
  gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 1 - gwas_geno1$maf_gwas, gwas_geno1$maf_gwas)
  gwas_geno1 <- gwas_geno1 %>%
    mutate(SNP = paste(chrbp, Base_unified, sep = ":"))
  
  # 合并GWAS数据
  coloc_data <- ASTHMA_geno1 %>% inner_join(gwas_geno1, by = "SNP")
  
  # 筛选目标区域数据
  coloc_file2 <- coloc_data %>%
    filter(CHR == chr_value & BP > bp_value - 500000 & BP < bp_value + 500000)
  
  # 进行共定位分析
  res <- jc_coloc(gene_value, coloc_file2)
  result <- rbind(result, res)
  
  # 画区域图
  outcome_coloc <- coloc_file2
  outcome_fn <- outcome_coloc %>% select(RSID, P_ASTHMA)
  names(outcome_fn) <- c('rsid', 'pval')
  gwas_fn <- outcome_coloc %>% select(RSID, P)
  names(gwas_fn) <- c('rsid', 'pval')
  
  # 找最显著的点
  data <- cbind(outcome_fn, gwas_fn[, 2])
  names(data) <- c('RSID', 'P_ASTHMA', 'P_phe')
  data$sum_p <- -log(data$P_ASTHMA) + -log(data$P_phe)
  rsid_max_sum_p <- data[which.max(data$sum_p), "RSID"]
  result2<-rbind(result2,rsid_max_sum_p)
  
  locus_plot <- paste0(gene_value, "_ASTHMA_", gsub("_assoc.result", "", filename), "_coloc.pdf")
  
  # 画图
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = outcome_fn, 
                    population = "EAS", 
                    title1 = "Phenotype", title2 = "ASTHMA", 
                    legend_position = 'topleft', 
                    marker_col1 = marker_col, pval_col1 = pval_col, marker_col2 = marker_col, pval_col2 = pval_col,
                    snp = rsid_max_sum_p,
                    genome = "hg19")+ theme(
    text = element_text(family = "Arial"),  
    plot.title = element_text(size = 15, family = "Arial"),  
    axis.title = element_text(size = 12, family = "Arial"), 
    axis.text = element_text(size = 10, family = "Arial"),  
    legend.text = element_text(size = 10, family = "Arial"), 
    legend.title = element_text(size = 12, family = "Arial") 
  )
  
  # 保存图像
  ggsave(locus_plot, plot = p, width = 9, height = 5, device = cairo_pdf)
  
}

# 输出最终结果
print(result)
print(result2)
result<-cbind(result,result2)
result$gene<-rownames(result)
result<-result[,c(8,1:7)]

write_xlsx(result, "/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ASTHMA/data/result_76loci_500kb_coloc.xlsx")






############################################ILD############################################
# define function
jc_coloc <- function(gene_id, coloc_file) {
  outcome_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_ILD,
    varbeta = coloc_file$varbeta_ILD,
    type = "cc",
    N = coloc_file$N,
    MAF = coloc_file$maf_ILD
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

result<-data.frame() #储存共定位结果
result2<-data.frame() #储存最显著的snp

data1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_sig_snps.clumped')[,c("CHR","BP","SNP")]
ILD_geno1<-fread('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/ILD_META_for_coloc.txt')
min_p_bolt_lmm_data_ILD<-subset(min_p_bolt_lmm_data,min_p_bolt_lmm_data$RSID %in% data1$SNP)

write.table(min_p_bolt_lmm_data_ILD,"/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/data/ILD_MATCH_phenotype.txt",quote=F,sep="\t",row.names=F)

# 设置图形保存路径
setwd('/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/data')
# 遍历 min_p_bolt_lmm_data_ILD 的1:10行
for (i in 1:nrow(min_p_bolt_lmm_data_ILD)) {
  
  # 获取当前行的信息
  current_row <- min_p_bolt_lmm_data_ILD[i, ]
  filename <- current_row$filename
  chr_value <- current_row$CHR
  bp_value <- current_row$BP
  gene_value <- current_row$gene 
  
  # 生成文件路径
  file_path <- file.path("/data1/Imaging_Assoc_2023/Assoc_Lung_phenotype_low_3.8w_clean", filename)
  
  # 读取文件
  gwas_geno <- fread(file_path)
  gwas_geno$chrbp <- paste(gwas_geno$CHR, gwas_geno$BP, sep = ":")
  
  # 重命名列
  names(gwas_geno)[names(gwas_geno) == 'SNP'] <- 'GWASsnpid'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE1'] <- 'A1'
  names(gwas_geno)[names(gwas_geno) == 'ALLELE0'] <- 'A2'
  names(gwas_geno)[names(gwas_geno) == 'A1FREQ'] <- 'maf_gwas'
  names(gwas_geno)[names(gwas_geno) == 'BETA'] <- 'beta_gwas'
  names(gwas_geno)[names(gwas_geno) == 'P_BOLT_LMM'] <- 'P'
  
  # 过滤和转换数据
  gwas_geno1 <- gwas_geno %>%
    filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
    mutate(Base_pair = paste(A1, A2, sep = ":")) %>%
    mutate(varbeta_gwas = (SE)^2) %>%
    select(CHR, BP, GWASsnpid, chrbp, Base_pair, maf_gwas, beta_gwas, varbeta_gwas, P, RSID)
  
  # 整理等位基因
  gwas_geno1$Base_unified <- ifelse(gwas_geno1$Base_pair %in% c("T:G", "G:T", "C:A"), "A:C", 
                                    ifelse(gwas_geno1$Base_pair %in% c("T:C", "C:T", "G:A"), "A:G",
                                    ifelse(gwas_geno1$Base_pair == "T:A", "A:T",
                                    ifelse(gwas_geno1$Base_pair == "G:C", "C:G", gwas_geno1$Base_pair))))
  gwas_geno1$beta_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), -gwas_geno1$beta_gwas, gwas_geno1$beta_gwas)
  gwas_geno1$maf_gwas <- ifelse(gwas_geno1$Base_pair %in% c("G:T", "C:A", "C:T", "G:A", "T:A", "G:C"), 1 - gwas_geno1$maf_gwas, gwas_geno1$maf_gwas)
  gwas_geno1 <- gwas_geno1 %>%
    mutate(SNP = paste(chrbp, Base_unified, sep = ":"))
  
  # 合并GWAS数据
  coloc_data <- ILD_geno1 %>% inner_join(gwas_geno1, by = "SNP")
  
  # 筛选目标区域数据
  coloc_file2 <- coloc_data %>%
    filter(CHR == chr_value & BP > bp_value - 500000 & BP < bp_value + 500000)
  
  # 进行共定位分析
  res <- jc_coloc(gene_value, coloc_file2)
  result <- rbind(result, res)
  
  # 画区域图
  outcome_coloc <- coloc_file2
  outcome_fn <- outcome_coloc %>% select(RSID, P_ILD)
  names(outcome_fn) <- c('rsid', 'pval')
  gwas_fn <- outcome_coloc %>% select(RSID, P)
  names(gwas_fn) <- c('rsid', 'pval')
  
  # 找最显著的点
  data <- cbind(outcome_fn, gwas_fn[, 2])
  names(data) <- c('RSID', 'P_ILD', 'P_phe')
  data$sum_p <- -log(data$P_ILD) + -log(data$P_phe)
  rsid_max_sum_p <- data[which.max(data$sum_p), "RSID"]
  result2<-rbind(result2,rsid_max_sum_p)
  
  locus_plot <- paste0(gene_value, "_ILD_", gsub("_assoc.result", "", filename), "_coloc.pdf")
  
  # 画图
  p <- locuscompare(in_fn1 = gwas_fn, in_fn2 = outcome_fn, 
                    population = "EAS", 
                    title1 = "Phenotype", title2 = "ILD", 
                    legend_position = 'topleft', 
                    marker_col1 = marker_col, pval_col1 = pval_col, marker_col2 = marker_col, pval_col2 = pval_col,
                    snp = rsid_max_sum_p,
                    genome = "hg19")+ theme(
    text = element_text(family = "Arial"),  
    plot.title = element_text(size = 15, family = "Arial"),  
    axis.title = element_text(size = 12, family = "Arial"), 
    axis.text = element_text(size = 10, family = "Arial"),  
    legend.text = element_text(size = 10, family = "Arial"), 
    legend.title = element_text(size = 12, family = "Arial") 
  )
  
  # 保存图像
  ggsave(locus_plot, plot = p, width = 9, height = 5, device = cairo_pdf)
  
}

# 输出最终结果
print(result)
print(result2)
result<-cbind(result,result2)
result$gene<-rownames(result)
result<-result[,c(8,1:7)]

write_xlsx(result, "/data1/Imaging_Assoc_2023/calculation_result_JC_new/7.coloc_JC/ILD/data/result_2loci_500kb_coloc.xlsx")

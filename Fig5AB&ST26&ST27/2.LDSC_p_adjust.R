library(openxlsx)

####
NSCLC <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 1)
COPD <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 2)
ILD <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 3)
ASTHMA <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 4)
RATIO <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 5)
FEV1 <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 6)
FVC <- read.xlsx('LDSC_RAW_res.xlsx',sheet = 7)

NSCLC$p.adj <- p.adjust(NSCLC$p,method ="fdr" )
COPD$p.adj <- p.adjust(COPD$p,method ="fdr" )
ILD$p.adj <- p.adjust(ILD$p,method ="fdr" )
ASTHMA$p.adj <- p.adjust(ASTHMA$p,method ="fdr" )
RATIO$p.adj <- p.adjust(RATIO$p,method ="fdr" )
FEV1$p.adj <- p.adjust(FEV1$p,method ="fdr" )
FVC$p.adj <- p.adjust(FVC$p,method ="fdr" )


write.table(NSCLC,"NSCLC_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(COPD,"COPD_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(ILD,"ILD_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(ASTHMA,"ASTHMA_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(RATIO,"RATIO_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(FEV1,"FEV1_adjusted.xlsx",quote = F,sep = "\t",row.names = F)
write.table(FVC,"FVC_adjusted.xlsx",quote = F,sep = "\t",row.names = F)

########################################################################################################################################
################################################查看显著位点在研究结局里的效应大小和方向#####################################################
########################################################################################################################################


R
library(data.table)

# 获取所有要提取的表型名称
test<-fread('160phe_unique_leadsnps-v2.txt') #all lead snps
test$CHR_BP<-paste(test$CHR,test$BP,sep=":")
outcome<-fread('NSCLC_for_merge_effect.txt') #同一bp位置的点在表型里的效应
colnames(outcome) <- paste0(colnames(outcome), "_LC")
outcome$CHR_BP<-paste(outcome$CHR_LC,outcome$BP_LC,sep=":")
outcome_2<-subset(outcome,outcome$CHR_BP %in% test$CHR_BP )

# 对CHR_BP列进行去重操作
outcome_2_unique <- outcome_2[!duplicated(outcome_2$CHR_BP), ]

test2<-merge(test,outcome_2_unique,by="CHR_BP")

#添加effect allele是否一致的一列(test代码没有),在这里NSCLC的beta是A1的效应，影像组summary数据的beta也是A1的效应,所以只需要比较两个A1的效应是否一致，一致就是原来的beta，不一致（翻转）则为-beta
test2$BETA_LC_new<- ifelse(test2$EA==test2$A1_LC & test2$NEA==test2$A2_LC, test2$BETA_LC, ifelse(test2$EA == test2$A2_LC & test2$NEA ==test2$A1_LC,-test2$BETA_LC,"Wrongsnp"))

# 将合并的数据保存到新的excel文件中
write.table(test2, "lead_snp_trait_clumping_merge_outcome.xlsx",, quote = FALSE, sep = '\t', row.names = FALSE )
#画MAPKAPK5的lead snp的区域图
#lead snp：12:112437149(rs11066122)

awk -v pos=112437149 -v flank=300000 'BEGIN{FS=","; OFS="\t"; print "RSID", "P_BOLT_LMM"} $2 == "12" && $3 >= (pos - flank) && $3 <= (pos + flank) {print $(NF), $16}' ./Assoc_Lung_phenotype_low_3.8w_clean/original_shape_Flatness_LL_assoc.result > "rs11066122.txt"
#plot in locus zoom


library("ggplot2")
library("ggbeeswarm")
library("ggsci")

rm(list=ls())

#加载协变量
cov=as.data.frame(t(read.table("338LungCombine.batch.combined_covariates_5pca.txt",h=T,check.names = F)))
cov$IID=row.names(cov)
names(cov)[1:54]=cov[1,c(1:54)]
cov=cov[-1,]
#转换数据类型
cov[, 1:54] <- lapply(cov[, 1:54], as.numeric)
cov$sex <- as.factor(cov$sex)
cov$smoking_status <- as.factor(cov$smoking_status)
cov$batch <- as.factor(cov$batch)


######################################################MAPKAPK5，rs11066122##########################################
MAPKAPK5=as.data.frame(t(read.table("338exp_for_MAPKAPK5.txt",h=T,check.names = F)))
MAPKAPK5$sampleid=row.names(MAPKAPK5)
names(MAPKAPK5)[1]=c("MAPKAPK5")
MAPKAPK5=MAPKAPK5[-1,]


geno=read.table("338_rs11066122.raw",h=T)
data=merge(geno,MAPKAPK5,by.x=c("IID"),by.y=c("sampleid"))
data$MAPKAPK5=as.numeric(data$MAPKAPK5)

data=merge(data,cov,by.x=c("IID"),by.y=c("IID"))

model <- glm(MAPKAPK5~pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

model2 <- glm(MAPKAPK5~X12.112437149.A.G_C+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

summary(model2) #1.90e-05

data$predict_MAPKAPK5 <- predict(model)
data$rs11066122 <- ifelse(data$X12.112437149.A.G_C==0,"TT",
                         ifelse(data$X12.112437149.A.G_C==1,"CT","CC"))

data$MAPKAPK5_plot <- data$MAPKAPK5-data$predict_MAPKAPK5

b4 = ggplot(data, aes(y = MAPKAPK5_plot, x = rs11066122, fill=rs11066122)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(color= rs11066122), alpha = 0.25, width = 0.15, color = "black",position = position_dodge(width = 0.5)) + 
  geom_violin(aes(x = rs11066122, fill = rs11066122), color = "lightgray", alpha = 0.5 , width = 0.5) +
  geom_beeswarm(aes(x = rs11066122, y = MAPKAPK5_plot, color = rs11066122, fill=rs11066122), alpha = 0.5, dodge.width=0.5)+
  scale_y_continuous(limits = c(-1.2,1.5), expand = c(0, 0)) +
  scale_fill_npg() +  scale_color_npg() + labs(x = "rs11066122",y = "MAPKAPK5 expression",fill = "")  + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour=NA), # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour=NA),  # 设置透明背景
        axis.line = element_line(colour = "black")
  ) +
  stat_summary(fun.y = "median", geom = "point", size = 1) + 
  stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1),color="gray")+
  annotate(geom = "text", x = 2, y = 1.4, label = "p = 1.90e-05", size = 4)

b4
ggsave("rs11066122_MAPKAPK5_eqtl.pdf", plot = b4, device = "pdf",width = 5,height = 5)


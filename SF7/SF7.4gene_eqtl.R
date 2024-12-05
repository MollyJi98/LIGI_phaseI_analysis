###################四个基因的区域图###################
library("ggplot2")
library("ggbeeswarm")
library("ggsci")

#QTL确认#
rm(list=ls())

cov=as.data.frame(t(read.table("338LungCombine.batch.combined_covariates_5pca.txt",h=T,check.names = F)))
cov$IID=row.names(cov)
names(cov)[1:54]=cov[1,c(1:54)]
cov=cov[-1,]
#转换数据类型
cov[, 1:54] <- lapply(cov[, 1:54], as.numeric)
cov$sex <- as.factor(cov$sex)
cov$smoking_status <- as.factor(cov$smoking_status)
cov$batch <- as.factor(cov$batch)

######################################################FAM13A，rs4505789##########################################
FAM13A=as.data.frame(t(read.table("338exp_for_FAM13A_HERC3.txt",h=T,check.names = F)[2,]))
FAM13A$sampleid=row.names(FAM13A)
names(FAM13A)[1]=c("FAM13A")
FAM13A=FAM13A[-1,]


geno=read.table("rs4505789.raw",h=T)
data=merge(geno,FAM13A,by.x=c("IID"),by.y=c("sampleid"))
data$FAM13A=as.numeric(data$FAM13A)

data=merge(data,cov,by.x=c("IID"),by.y=c("IID"))


#画图#

model <- glm(FAM13A~pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

model2 <- glm(FAM13A~X4.89825615.A.C_A+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

summary(model2) #2.90e-14

data$predict_FAM13A <- predict(model)
data$rs4505789 <- ifelse(data$X4.89825615.A.C_A==0,"CC",
                         ifelse(data$X4.89825615.A.C_A==1,"CA","AA"))

data$FAM13A_plot <- data$FAM13A-data$predict_FAM13A



b4 = ggplot(data, aes(y = FAM13A_plot, x = rs4505789, fill=rs4505789)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(color= rs4505789), alpha = 0.25, width = 0.15, color = "black",position = position_dodge(width = 0.5)) + 
  geom_violin(aes(x = rs4505789, fill = rs4505789), color = "lightgray", alpha = 0.5 , width = 0.5) +
  geom_beeswarm(aes(x = rs4505789, y = FAM13A_plot, color = rs4505789, fill=rs4505789), alpha = 0.5, dodge.width=0.5)+
  scale_y_continuous(limits = c(-1.2,1.5), expand = c(0, 0)) +
  scale_fill_npg() +  scale_color_npg() + labs(x = "rs4505789",y = "FAM13A expression",fill = "")  + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour=NA), # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour=NA),  # 设置透明背景
        axis.line = element_line(colour = "black")
  ) +
  stat_summary(fun.y = "median", geom = "point", size = 1) + 
  stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1),color="gray")+
  annotate(geom = "text", x = 2, y = 1.2, label = "p = 2.90e-14", size = 4)

b4
ggsave("rs4505789_FAM13A_eqtl.pdf", plot = b4, device = "pdf",width = 5,height = 5)




######################################################AGER,s10947233##########################################
AGER=as.data.frame(t(read.table("338exp_for_AGER.txt",h=T,check.names = F)))
AGER$sampleid=row.names(AGER)
names(AGER)[1]=c("AGER")
AGER=AGER[-1,]


geno=read.table("338_rs10947233.raw",h=T)
data=merge(geno,AGER,by.x=c("IID"),by.y=c("sampleid"))
data$AGER=as.numeric(data$AGER)

data=merge(data,cov,by.x=c("IID"),by.y=c("IID"))


#画图#

model <- glm(AGER~pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

model2 <- glm(AGER~X6.32124424.A.C_T+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)
summary(model2) #2.81e-10

data$predict_AGER <- predict(model)
data$rs10947233 <- ifelse(data$X6.32124424.A.C_T==0,"GG",
                         ifelse(data$X6.32124424.A.C_T==1,"GT","TT"))

data$AGER_plot <- data$AGER-data$predict_AGER


b4 = ggplot(data, aes(y = AGER_plot, x = rs10947233, fill=rs10947233)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(color= rs10947233), alpha = 0.25, width = 0.15, color = "black",position = position_dodge(width = 0.5)) + 
  geom_violin(aes(x = rs10947233, fill = rs10947233), color = "lightgray", alpha = 0.5 , width = 0.5) +
  geom_beeswarm(aes(x = rs10947233, y = AGER_plot, color = rs10947233, fill=rs10947233), alpha = 0.5, dodge.width=0.5)+
  scale_y_continuous(limits = c(-1.5,1.2), expand = c(0, 0)) +
  scale_fill_npg() +  scale_color_npg() + labs(x = "rs10947233",y = "AGER expression",fill = "")  + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour=NA), # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour=NA),  # 设置透明背景
        axis.line = element_line(colour = "black")
  ) +
  stat_summary(fun.y = "median", geom = "point", size = 1) + 
  stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1),color="gray")+
  annotate(geom = "text", x = 2, y = 1.1, label = "p = 2.81e-10", size = 4)

b4
ggsave("rs10947233_AGER_eqtl.pdf", plot = b4, device = "pdf",width = 5,height = 5)





######################################################MFAP2 rs6682671##########################################
MFAP2=as.data.frame(t(read.table("338exp_for_MFAP2.txt",h=T,check.names = F)))
MFAP2$sampleid=row.names(MFAP2)
names(MFAP2)[1]=c("MFAP2")
MFAP2=MFAP2[-1,]


geno=read.table("338_rs6682671.raw",h=T)
data=merge(geno,MFAP2,by.x=c("IID"),by.y=c("sampleid"))
data$MFAP2=as.numeric(data$MFAP2)

data=merge(data,cov,by.x=c("IID"),by.y=c("IID"))


#画图#

model <- glm(MFAP2~pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

model2 <- glm(MFAP2~X1.17311882.A.G_T+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)
summary(model2) #5.10e-10

data$predict_MFAP2 <- predict(model)
data$rs6682671 <- ifelse(data$X1.17311882.A.G_T==0,"CC",
                         ifelse(data$X1.17311882.A.G_T==1,"CT","TT"))

data$MFAP2_plot <- data$MFAP2-data$predict_MFAP2
summary(data$MFAP2_plot )

b4 = ggplot(data, aes(y = MFAP2_plot, x = rs6682671, fill=rs6682671)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(color= rs6682671), alpha = 0.25, width = 0.15, color = "black",position = position_dodge(width = 0.5)) + 
  geom_violin(aes(x = rs6682671, fill = rs6682671), color = "lightgray", alpha = 0.5 , width = 0.5) +
  geom_beeswarm(aes(x = rs6682671, y = MFAP2_plot, color = rs6682671, fill=rs6682671), alpha = 0.5, dodge.width=0.5)+
  scale_y_continuous(limits = c(-1.7,1.6), expand = c(0, 0)) +
  scale_fill_npg() +  scale_color_npg() + labs(x = "rs6682671",y = "MFAP2 expression",fill = "")  + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour=NA), # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour=NA),  # 设置透明背景
        axis.line = element_line(colour = "black")
  ) +
  stat_summary(fun.y = "median", geom = "point", size = 1) + 
  stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1),color="gray")+
  annotate(geom = "text", x = 2, y = 1.4, label = "p = 5.10e-10", size = 4)

b4
ggsave("rs6682671_MFAP2_eqtl.pdf", plot = b4, device = "pdf",width = 5,height = 5)




######################################################ADGRG6（GPR126）rs6570508##########################################
GPR126=as.data.frame(t(read.table("338exp_for_GPR126.txt",h=T,check.names = F)))
GPR126$sampleid=row.names(GPR126)
names(GPR126)[1]=c("GPR126")
GPR126=GPR126[-1,]


geno=read.table("338_rs6570508.raw",h=T)
data=merge(geno,GPR126,by.x=c("IID"),by.y=c("sampleid"))
data$GPR126=as.numeric(data$GPR126)

data=merge(data,cov,by.x=c("IID"),by.y=c("IID"))


#画图#

model <- glm(GPR126~pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

model2 <- glm(GPR126~X6.142713842.A.G_G+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)
summary(model2) #0.031990
 

data$predict_GPR126 <- predict(model)
data$rs6570508 <- ifelse(data$X6.142713842.A.G_G==0,"AA",
                         ifelse(data$X6.142713842.A.G_G==1,"AG","GG"))

data$GPR126_plot <- data$GPR126-data$predict_GPR126
summary(data$GPR126_plot )

b4 = ggplot(data, aes(y = GPR126_plot, x = rs6570508, fill=rs6570508)) +
  stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(color= rs6570508), alpha = 0.25, width = 0.15, color = "black",position = position_dodge(width = 0.5)) + 
  geom_violin(aes(x = rs6570508, fill = rs6570508), color = "lightgray", alpha = 0.5 , width = 0.5) +
  geom_beeswarm(aes(x = rs6570508, y = GPR126_plot, color = rs6570508, fill=rs6570508), alpha = 0.5, dodge.width=0.5)+
  scale_y_continuous(limits = c(-1.3,0.9), expand = c(0, 0)) +
  scale_fill_npg() +  scale_color_npg() + labs(x = "rs6570508",y = "ADGRG6 expression",fill = "")  + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour=NA), # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour=NA),  # 设置透明背景
        axis.line = element_line(colour = "black")
  ) +
  stat_summary(fun.y = "median", geom = "point", size = 1) + 
  stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1),color="gray")+
  annotate(geom = "text", x = 2, y = 0.8, label = "p = 0.032", size = 4)

b4
ggsave("rs6570508_GPR126_eqtl.pdf", plot = b4, device = "pdf",width = 5,height = 5)

###################same for eqtl plot in GTEx data###################
rm(list=ls())

#调用R包
library(readxl)
library(ggthemes)
library(ggplot2)
library(ggpubr)

exp <- read_excel("准备文件_LQ.xlsx",sheet="表型并集")
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###TARGET_OUTCOME### ### ### ### ### ### ### ### ### ###
TARGET_OUTCOME <- read_excel("准备文件_LQ.xlsx",sheet="TARGET_OUTCOME")[,c(4:6,8:10)] #提取需要画图的数据
TARGET_OUTCOME_1 <- TARGET_OUTCOME[TARGET_OUTCOME$method=="Inverse variance weighted",]
TARGET_OUTCOME_2 <- merge(TARGET_OUTCOME_1,exp,all=T)
TARGET_OUTCOME_3 <- TARGET_OUTCOME_2[!is.na(TARGET_OUTCOME_2$check),]
TARGET_OUTCOME_3$OR <- round(exp(TARGET_OUTCOME_3$b),2)
TARGET_OUTCOME_3$lower_ci <- round(TARGET_OUTCOME_3$b-1.96*TARGET_OUTCOME_3$se,2)
TARGET_OUTCOME_3$upper_ci <- round(TARGET_OUTCOME_3$b+1.96*TARGET_OUTCOME_3$se,2)
TARGET_OUTCOME_3$pval_sig <- ifelse(TARGET_OUTCOME_3$pval<0.05,1,0)
TARGET_OUTCOME_3$padj_sig <- ifelse(TARGET_OUTCOME_3$pval<0.05/128,1,0)
###绘图
TARGET_OUTCOME_3$pval_sig <- as.factor(TARGET_OUTCOME_3$pval_sig)
TARGET_OUTCOME_3$n <- paste(TARGET_OUTCOME_3$exposure,TARGET_OUTCOME_3$Lobe,sep = "_")

desired_order <- c("10Percentile_CR","10Percentile_LL","10Percentile_LR","10Percentile_UL","10Percentile_UR",
                   "90Percentile_CR","90Percentile_LL","90Percentile_LR","90Percentile_UL","90Percentile_UR",
                   "Energy_LL","Entropy_CR","Entropy_LR","Entropy_UL","Entropy_UR",
                   "InterquartileRange_UL","InterquartileRange_UR","Kurtosis_CR","Kurtosis_LR","Kurtosis_UR",
                   "Mean_CR","Mean_LL","Mean_LR","Mean_UL","Mean_UR",
                   "MeanAbsoluteDeviation_UL","Median_CR","Median_LL","Median_LR","Median_UL","Median_UR",
                   "RobustMeanAbsoluteDeviation_UL","RobustMeanAbsoluteDeviation_UR","RootMeanSquared_CR","RootMeanSquared_LL","RootMeanSquared_LR",
                   "RootMeanSquared_UL","RootMeanSquared_UR","Skewness_CR","Skewness_LL","Skewness_LR",
                   "Skewness_UL","Skewness_UR","TotalEnergy_LL","Uniformity_CR","Uniformity_LR",
                   "Uniformity_UL","Uniformity_UR","Variance_UR","Elongation_CR",
                   "LeastAxisLength_CR","LeastAxisLength_UL","LeastAxisLength_UR","MajorAxisLength_CR","MajorAxisLength_UL",
                   "Maximum2DDiameterColumn_UR","Maximum2DDiameterRow_UR","Maximum3DDiameter_LR","Maximum3DDiameter_UL","Maximum3DDiameter_UR",
                   "MeshVolume_LL","MeshVolume_UL","MeshVolume_UR","MinorAxisLength_LL","MinorAxisLength_UR",
                   "Sphericity_UR","SurfaceArea_CR","SurfaceArea_UL","SurfaceArea_UR",
                   "SurfaceVolumeRatio_LR","VoxelVolume_LL","VoxelVolume_UL","VoxelVolume_UR")
TARGET_OUTCOME_3 <- TARGET_OUTCOME_3[order(match(TARGET_OUTCOME_3$n, desired_order)), ]
TARGET_OUTCOME_3$n <- factor(TARGET_OUTCOME_3$n, levels = rev(desired_order))
max(TARGET_OUTCOME_3$upper_ci) #3.05
min(TARGET_OUTCOME_3$lower_ci) #-3.65
p5<-ggplot(TARGET_OUTCOME_3,aes(color = pval_sig,x = b , xmin = lower_ci, xmax = upper_ci , y = n))+
  geom_pointrange(size=0,fatten = 0)+ 
  scale_color_manual(values = c("#A0A0A4", "#A81D2E"))+ #设置颜色
  theme_bw()+ #设置经典主题
  coord_cartesian(xlim = c(-4.0, 4.0))+ #控制显示范围
  xlab("")+ #X轴标题
  ylab("")+ #Y轴标题，此处设置为空值
  ggtitle("TARGET_OUTCOME")+ #主标题
  theme(
    legend.position = "none", 
    axis.text.y = element_text(color = "black"), #设置Y轴字体颜色
    axis.ticks.x = element_blank(), # 删除X轴刻度
    axis.ticks.y = element_blank(), # 删除Y轴刻度
    plot.title = element_text(size = 10 , face = "bold") #设置标题大小
  )+
  #添加背景颜色
  annotate("rect", xmin = -4.5, xmax = 4.5, ymin = 0, ymax = 24, alpha = .5, fill = "#fef4f3")+ #添加背景色块
  annotate("rect", xmin = -4.5, xmax = 4.5, ymin = 24, ymax = 74, alpha = .5, fill = "#e0edf9")+ #添加背景色块
  scale_x_continuous(breaks = c(-3.0,-1.5,0,1.5,3.0))+ #设置X轴刻度
  geom_pointrange(
    linewidth=0.5, # 控制线的宽度
    fatten = 0.6)+ 
  geom_vline(xintercept = 0, color = "black",linewidth = 0.4)+ #添加直线
  # 添加 * 号：变量 padj_sig 为 1 的位置
  geom_text(
    data = TARGET_OUTCOME_3[TARGET_OUTCOME_3$padj_sig == 1, ], 
    aes(x = upper_ci + 0.1, y = n, label = "*"), 
    color = "black", 
    size = 4)
#保存图像
ggsave("TARGET_OUTCOME.pdf",plot=p5,width=6,height=12)



###################画完后在ai里拼图###################




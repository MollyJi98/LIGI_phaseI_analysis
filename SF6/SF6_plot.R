rm(list=ls())

library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

###########################################first order###########################################
f1 <- read_excel("准备文件.xlsx",sheet="first order") #文件是整理好的每个遗传区域在不同表型里的显著情况，整理的结果来源于ST5,7的显著表型数量
f_varname <- c('10Percentile_CR','10Percentile_UR','10Percentile_UL','10Percentile_LR','10Percentile_LL','90Percentile_CR','90Percentile_UR','90Percentile_UL','90Percentile_LR','90Percentile_LL','Energy_CR','Energy_UR','Energy_UL','Energy_LR','Energy_LL','Entropy_CR','Entropy_UR','Entropy_UL','Entropy_LR','Entropy_LL','InterquartileRange_CR','InterquartileRange_UR','InterquartileRange_UL','InterquartileRange_LR','InterquartileRange_LL','Kurtosis_CR','Kurtosis_UR','Kurtosis_UL','Kurtosis_LR','Kurtosis_LL','Maximum_CR','Maximum_UR','Maximum_UL','Maximum_LR','Maximum_LL','Mean_CR','Mean_UR','Mean_UL','Mean_LR','Mean_LL','MeanAbsoluteDeviation_CR','MeanAbsoluteDeviation_UR','MeanAbsoluteDeviation_UL','MeanAbsoluteDeviation_LR','MeanAbsoluteDeviation_LL','Median_CR','Median_UR','Median_UL','Median_LR','Median_LL','Minimum_CR','Minimum_UR','Minimum_UL','Minimum_LR','Minimum_LL','Range_CR','Range_UR','Range_UL','Range_LR','Range_LL','RobustMeanAbsoluteDeviation_CR','RobustMeanAbsoluteDeviation_UR','RobustMeanAbsoluteDeviation_UL','RobustMeanAbsoluteDeviation_LR','RobustMeanAbsoluteDeviation_LL','RootMeanSquared_CR','RootMeanSquared_UR','RootMeanSquared_UL','RootMeanSquared_LR','RootMeanSquared_LL','Skewness_CR','Skewness_UR','Skewness_UL','Skewness_LR','Skewness_LL','TotalEnergy_CR','TotalEnergy_UR','TotalEnergy_UL','TotalEnergy_LR','TotalEnergy_LL','Uniformity_CR','Uniformity_UR','Uniformity_UL','Uniformity_LR','Uniformity_LL','Variance_CR','Variance_UR','Variance_UL','Variance_LR','Variance_LL')
# 循环遍历f_varname中的每个值
for (i in seq_along(f_varname)) {
  # 创建新变量的名称
  new_var <- f_varname[i]
  # 检查phenotype中是否包含当前值，并赋值给新变量
  f1[[new_var]] <- ifelse(grepl(paste0("\\b", f_varname[i], "\\b"), f1$phenotype), 1, 0)
}

f2 <- f1[,c(3,4,7:96)]

#####按照cytoBand排序
# 提取数字和字母部分
extracted1 <- str_match(f2$cytoBand, "(\\d+)([pq])(\\d+\\.?\\d*)")
# 创建数据框
extracted_f2 <- data.frame(
  cytoBand = f2$cytoBand,
  num1 = as.numeric(extracted1[, 2]),  # 第一个数字
  letter = extracted1[, 3],             # p或q
  num2 = as.numeric(extracted1[, 4])    # 第二个数字
)
# 按照排序规则排序
sorted_index <- order(extracted_f2$num1, extracted_f2$letter, extracted_f2$num2)
# 使用排序后的索引对数据框进行重新排序
f3 <- f2[sorted_index, ]
# 将数据变换成长格式
f3_long <- pivot_longer(f3, cols = c('10Percentile_CR','10Percentile_UR','10Percentile_UL','10Percentile_LR','10Percentile_LL','90Percentile_CR','90Percentile_UR','90Percentile_UL','90Percentile_LR','90Percentile_LL','Energy_CR','Energy_UR','Energy_UL','Energy_LR','Energy_LL','Entropy_CR','Entropy_UR','Entropy_UL','Entropy_LR','Entropy_LL','InterquartileRange_CR','InterquartileRange_UR','InterquartileRange_UL','InterquartileRange_LR','InterquartileRange_LL','Kurtosis_CR','Kurtosis_UR','Kurtosis_UL','Kurtosis_LR','Kurtosis_LL','Maximum_CR','Maximum_UR','Maximum_UL','Maximum_LR','Maximum_LL','Mean_CR','Mean_UR','Mean_UL','Mean_LR','Mean_LL','MeanAbsoluteDeviation_CR','MeanAbsoluteDeviation_UR','MeanAbsoluteDeviation_UL','MeanAbsoluteDeviation_LR','MeanAbsoluteDeviation_LL','Median_CR','Median_UR','Median_UL','Median_LR','Median_LL','Minimum_CR','Minimum_UR','Minimum_UL','Minimum_LR','Minimum_LL','Range_CR','Range_UR','Range_UL','Range_LR','Range_LL','RobustMeanAbsoluteDeviation_CR','RobustMeanAbsoluteDeviation_UR','RobustMeanAbsoluteDeviation_UL','RobustMeanAbsoluteDeviation_LR','RobustMeanAbsoluteDeviation_LL','RootMeanSquared_CR','RootMeanSquared_UR','RootMeanSquared_UL','RootMeanSquared_LR','RootMeanSquared_LL','Skewness_CR','Skewness_UR','Skewness_UL','Skewness_LR','Skewness_LL','TotalEnergy_CR','TotalEnergy_UR','TotalEnergy_UL','TotalEnergy_LR','TotalEnergy_LL','Uniformity_CR','Uniformity_UR','Uniformity_UL','Uniformity_LR','Uniformity_LL','Variance_CR','Variance_UR','Variance_UL','Variance_LR','Variance_LL'),
                        names_to = "phenotype", values_to = "Value")
#更改横轴和纵轴的排列顺序
f3_long$phenotype <- factor(f3_long$phenotype, levels = c('10Percentile_CR','10Percentile_UR','10Percentile_UL','10Percentile_LR','10Percentile_LL','90Percentile_CR','90Percentile_UR','90Percentile_UL','90Percentile_LR','90Percentile_LL','Energy_CR','Energy_UR','Energy_UL','Energy_LR','Energy_LL','Entropy_CR','Entropy_UR','Entropy_UL','Entropy_LR','Entropy_LL','InterquartileRange_CR','InterquartileRange_UR','InterquartileRange_UL','InterquartileRange_LR','InterquartileRange_LL','Kurtosis_CR','Kurtosis_UR','Kurtosis_UL','Kurtosis_LR','Kurtosis_LL','Maximum_CR','Maximum_UR','Maximum_UL','Maximum_LR','Maximum_LL','Mean_CR','Mean_UR','Mean_UL','Mean_LR','Mean_LL','MeanAbsoluteDeviation_CR','MeanAbsoluteDeviation_UR','MeanAbsoluteDeviation_UL','MeanAbsoluteDeviation_LR','MeanAbsoluteDeviation_LL','Median_CR','Median_UR','Median_UL','Median_LR','Median_LL','Minimum_CR','Minimum_UR','Minimum_UL','Minimum_LR','Minimum_LL','Range_CR','Range_UR','Range_UL','Range_LR','Range_LL','RobustMeanAbsoluteDeviation_CR','RobustMeanAbsoluteDeviation_UR','RobustMeanAbsoluteDeviation_UL','RobustMeanAbsoluteDeviation_LR','RobustMeanAbsoluteDeviation_LL','RootMeanSquared_CR','RootMeanSquared_UR','RootMeanSquared_UL','RootMeanSquared_LR','RootMeanSquared_LL','Skewness_CR','Skewness_UR','Skewness_UL','Skewness_LR','Skewness_LL','TotalEnergy_CR','TotalEnergy_UR','TotalEnergy_UL','TotalEnergy_LR','TotalEnergy_LL','Uniformity_CR','Uniformity_UR','Uniformity_UL','Uniformity_LR','Uniformity_LL','Variance_CR','Variance_UR','Variance_UL','Variance_LR','Variance_LL'))
f3_long$cytoBand_nearestgene <- factor(f3_long$cytoBand_nearestgene, levels = rev(unique(f3_long$cytoBand_nearestgene)))
# 创建热图
p1 <- ggplot(f3_long, aes(x = phenotype, y = cytoBand_nearestgene, fill = factor(Value))) +
  geom_tile(color = "#737373", width = 1, height = 1) +
  scale_fill_manual(values = c("white", "skyblue"), labels = c("0" = "White", "1" = "Blue")) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Value") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 5, vjust = 1, angle = 45,hjust=1), # 设置横轴标签字号
    axis.text.y = element_text(size = 5), # 设置纵轴标签字号
    axis.title = element_text(size = 10),  # 设置轴标题字号
    plot.title = element_text(size = 10),  # 设置图形标题字号
    panel.grid.major = element_blank(),     # 隐藏网格线
    panel.grid.minor = element_blank(),     # 隐藏网格线
    legend.position = "none"               # 不显示图例
  ) 
ggsave("plot_fistorder.pdf", plot = p1, width = 8, height = 3.9)







###########################################shape###########################################
s1 <- read_excel("准备文件.xlsx",sheet="shape")
s_varname <- c('Elongation_CR','Elongation_UR','Elongation_UL','Elongation_LR','Elongation_LL','Flatness_CR','Flatness_UR','Flatness_UL','Flatness_LR','Flatness_LL','LeastAxisLength_CR','LeastAxisLength_UR','LeastAxisLength_UL','LeastAxisLength_LR','LeastAxisLength_LL','MajorAxisLength_CR','MajorAxisLength_UR','MajorAxisLength_UL','MajorAxisLength_LR','MajorAxisLength_LL','Maximum2DDiameterColumn_CR','Maximum2DDiameterColumn_UR','Maximum2DDiameterColumn_UL','Maximum2DDiameterColumn_LR','Maximum2DDiameterColumn_LL','Maximum2DDiameterRow_CR','Maximum2DDiameterRow_UR','Maximum2DDiameterRow_UL','Maximum2DDiameterRow_LR','Maximum2DDiameterRow_LL','Maximum2DDiameterSlice_CR','Maximum2DDiameterSlice_UR','Maximum2DDiameterSlice_UL','Maximum2DDiameterSlice_LR','Maximum2DDiameterSlice_LL','Maximum3DDiameter_CR','Maximum3DDiameter_UR','Maximum3DDiameter_UL','Maximum3DDiameter_LR','Maximum3DDiameter_LL','MeshVolume_CR','MeshVolume_UR','MeshVolume_UL','MeshVolume_LR','MeshVolume_LL','MinorAxisLength_CR','MinorAxisLength_UR','MinorAxisLength_UL','MinorAxisLength_LR','MinorAxisLength_LL','Sphericity_CR','Sphericity_UR','Sphericity_UL','Sphericity_LR','Sphericity_LL','SurfaceArea_CR','SurfaceArea_UR','SurfaceArea_UL','SurfaceArea_LR','SurfaceArea_LL','SurfaceVolumeRatio_CR','SurfaceVolumeRatio_UR','SurfaceVolumeRatio_UL','SurfaceVolumeRatio_LR','SurfaceVolumeRatio_LL','VoxelVolume_CR','VoxelVolume_UR','VoxelVolume_UL','VoxelVolume_LR','VoxelVolume_LL')
# 循环遍历s_varname中的每个值
for (i in seq_along(s_varname)) {
  # 创建新变量的名称
  new_var <- s_varname[i]
  # 检查phenotype中是否包含当前值，并赋值给新变量
  s1[[new_var]] <- ifelse(grepl(paste0("\\b", s_varname[i], "\\b"), s1$phenotype), 1, 0)
}

s2 <- s1[,c(3,4,7:76)]

#####按照cytoBand排序
# 提取数字和字母部分
extracted2 <- str_match(s2$cytoBand, "(\\d+)([pq])(\\d+\\.?\\d*)")
# 创建数据框
extracted_s2 <- data.frame(
  cytoBand = s2$cytoBand,
  num1 = as.numeric(extracted2[, 2]),  # 第一个数字
  letter = extracted2[, 3],             # p或q
  num2 = as.numeric(extracted2[, 4])    # 第二个数字
)
# 按照排序规则排序
sorted_index2 <- order(extracted_s2$num1, extracted_s2$letter, extracted_s2$num2)
# 使用排序后的索引对数据框进行重新排序
s3 <- s2[sorted_index2, ]
# 将数据变换成长格式
s3_long <- pivot_longer(s3, cols = c('Elongation_CR','Elongation_UR','Elongation_UL','Elongation_LR','Elongation_LL','Flatness_CR','Flatness_UR','Flatness_UL','Flatness_LR','Flatness_LL','LeastAxisLength_CR','LeastAxisLength_UR','LeastAxisLength_UL','LeastAxisLength_LR','LeastAxisLength_LL','MajorAxisLength_CR','MajorAxisLength_UR','MajorAxisLength_UL','MajorAxisLength_LR','MajorAxisLength_LL','Maximum2DDiameterColumn_CR','Maximum2DDiameterColumn_UR','Maximum2DDiameterColumn_UL','Maximum2DDiameterColumn_LR','Maximum2DDiameterColumn_LL','Maximum2DDiameterRow_CR','Maximum2DDiameterRow_UR','Maximum2DDiameterRow_UL','Maximum2DDiameterRow_LR','Maximum2DDiameterRow_LL','Maximum2DDiameterSlice_CR','Maximum2DDiameterSlice_UR','Maximum2DDiameterSlice_UL','Maximum2DDiameterSlice_LR','Maximum2DDiameterSlice_LL','Maximum3DDiameter_CR','Maximum3DDiameter_UR','Maximum3DDiameter_UL','Maximum3DDiameter_LR','Maximum3DDiameter_LL','MeshVolume_CR','MeshVolume_UR','MeshVolume_UL','MeshVolume_LR','MeshVolume_LL','MinorAxisLength_CR','MinorAxisLength_UR','MinorAxisLength_UL','MinorAxisLength_LR','MinorAxisLength_LL','Sphericity_CR','Sphericity_UR','Sphericity_UL','Sphericity_LR','Sphericity_LL','SurfaceArea_CR','SurfaceArea_UR','SurfaceArea_UL','SurfaceArea_LR','SurfaceArea_LL','SurfaceVolumeRatio_CR','SurfaceVolumeRatio_UR','SurfaceVolumeRatio_UL','SurfaceVolumeRatio_LR','SurfaceVolumeRatio_LL','VoxelVolume_CR','VoxelVolume_UR','VoxelVolume_UL','VoxelVolume_LR','VoxelVolume_LL'),
                        names_to = "phenotype", values_to = "Value")
#更改横轴和纵轴的排列顺序
s3_long$phenotype <- factor(s3_long$phenotype, levels = c('Elongation_CR','Elongation_UR','Elongation_UL','Elongation_LR','Elongation_LL','Flatness_CR','Flatness_UR','Flatness_UL','Flatness_LR','Flatness_LL','LeastAxisLength_CR','LeastAxisLength_UR','LeastAxisLength_UL','LeastAxisLength_LR','LeastAxisLength_LL','MajorAxisLength_CR','MajorAxisLength_UR','MajorAxisLength_UL','MajorAxisLength_LR','MajorAxisLength_LL','Maximum2DDiameterColumn_CR','Maximum2DDiameterColumn_UR','Maximum2DDiameterColumn_UL','Maximum2DDiameterColumn_LR','Maximum2DDiameterColumn_LL','Maximum2DDiameterRow_CR','Maximum2DDiameterRow_UR','Maximum2DDiameterRow_UL','Maximum2DDiameterRow_LR','Maximum2DDiameterRow_LL','Maximum2DDiameterSlice_CR','Maximum2DDiameterSlice_UR','Maximum2DDiameterSlice_UL','Maximum2DDiameterSlice_LR','Maximum2DDiameterSlice_LL','Maximum3DDiameter_CR','Maximum3DDiameter_UR','Maximum3DDiameter_UL','Maximum3DDiameter_LR','Maximum3DDiameter_LL','MeshVolume_CR','MeshVolume_UR','MeshVolume_UL','MeshVolume_LR','MeshVolume_LL','MinorAxisLength_CR','MinorAxisLength_UR','MinorAxisLength_UL','MinorAxisLength_LR','MinorAxisLength_LL','Sphericity_CR','Sphericity_UR','Sphericity_UL','Sphericity_LR','Sphericity_LL','SurfaceArea_CR','SurfaceArea_UR','SurfaceArea_UL','SurfaceArea_LR','SurfaceArea_LL','SurfaceVolumeRatio_CR','SurfaceVolumeRatio_UR','SurfaceVolumeRatio_UL','SurfaceVolumeRatio_LR','SurfaceVolumeRatio_LL','VoxelVolume_CR','VoxelVolume_UR','VoxelVolume_UL','VoxelVolume_LR','VoxelVolume_LL'))
s3_long$cytoBand_nearestgene <- factor(s3_long$cytoBand_nearestgene, levels = rev(unique(s3_long$cytoBand_nearestgene)))
# 创建热图
p2 <- ggplot(s3_long, aes(x = phenotype, y = cytoBand_nearestgene, fill = factor(Value))) +
  geom_tile(color = "#737373", width = 1, height = 1) +
  scale_fill_manual(values = c("white", "skyblue"), labels = c("0" = "White", "1" = "Blue")) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Value") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 5, vjust = 1, angle = 45,hjust=1), # 设置横轴标签字号
    axis.text.y = element_text(size = 5), # 设置纵轴标签字号
    axis.title = element_text(size = 10),  # 设置轴标题字号
    plot.title = element_text(size = 10),  # 设置图形标题字号
    panel.grid.major = element_blank(),     # 隐藏网格线
    panel.grid.minor = element_blank(),     # 隐藏网格线
    legend.position = "none"               # 不显示图例
  )

ggsave("plot_shape.pdf", plot = p2, width = 8, height = 11.6)


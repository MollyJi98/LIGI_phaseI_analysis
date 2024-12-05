########基于LDSC的结果计算截距项########
library(ggplot2)
library(gridExtra)

# 定义文件路径
file_path <- "./step3_calculate_LDSC/ILD_meta.log"

# 读取文件内容
file_content <- readLines(file_path)

# 初始化一个空向量来存储结果
intercept_content <- c()

# 提取每个 "Intercept:" 行后面直到括号 "(" 之前的内容
for (line in file_content) {
  if (grepl("Intercept:", line)) {
    intercept <- sub(".*Intercept: (.*?)\\(.*", "\\1", line)
    intercept_content <- c(intercept_content, intercept)
  }
}

intercept_content <- as.data.frame(intercept_content)
# 删除第一行
intercept_content <- intercept_content[-1, , drop = FALSE]

# 只保留奇数行
intercept_content <- intercept_content[seq(1, nrow(intercept_content), by = 2), , drop = FALSE]

#colnames定义
names(intercept_content)[1] <-'Intercept'
# names <- c("10Percentile_CR","10Percentile_LL","10Percentile_LR","10Percentile_UL","10Percentile_UR","90Percentile_CR","90Percentile_LL","90Percentile_LR","90Percentile_UL","90Percentile_UR","Energy_CR","Energy_LL","Energy_LR","Energy_UL","Energy_UR","Entropy_CR","Entropy_LL","Entropy_LR","Entropy_UL","Entropy_UR","InterquartileRange_CR","InterquartileRange_LL","InterquartileRange_LR","InterquartileRange_UL","InterquartileRange_UR","Kurtosis_CR","Kurtosis_LL","Kurtosis_LR","Kurtosis_UL","Kurtosis_UR","Maximum_CR","Maximum_LL","Maximum_LR","Maximum_UL","Maximum_UR","Mean_CR","Mean_LL","Mean_LR","Mean_UL","Mean_UR","MeanAbsoluteDeviation_CR","MeanAbsoluteDeviation_LL","MeanAbsoluteDeviation_LR","MeanAbsoluteDeviation_UL","MeanAbsoluteDeviation_UR","Median_CR","Median_LL","Median_LR","Median_UL","Median_UR","Minimum_CR","Minimum_LL","Minimum_LR","Minimum_UL","Range_CR","Range_LL","Range_LR","Range_UL","Range_UR","RobustMeanAbsoluteDeviation_CR","RobustMeanAbsoluteDeviation_LL","RobustMeanAbsoluteDeviation_LR","RobustMeanAbsoluteDeviation_UL","RobustMeanAbsoluteDeviation_UR","RootMeanSquared_CR","RootMeanSquared_LL","RootMeanSquared_LR","RootMeanSquared_UL","RootMeanSquared_UR","Skewness_CR","Skewness_LL","Skewness_LR","Skewness_UL","Skewness_UR","TotalEnergy_CR","TotalEnergy_LL","TotalEnergy_LR","TotalEnergy_UL","TotalEnergy_UR","Uniformity_CR","Uniformity_LL","Uniformity_LR","Uniformity_UL","Uniformity_UR","Variance_CR","Variance_LL","Variance_LR","Variance_UL","Variance_UR","Elongation_CR","Elongation_LL","Elongation_LR","Elongation_UL","Elongation_UR","Flatness_CR","Flatness_LL","Flatness_LR","Flatness_UL","Flatness_UR","LeastAxisLength_CR","LeastAxisLength_LL","LeastAxisLength_LR","LeastAxisLength_UL","LeastAxisLength_UR","MajorAxisLength_CR","MajorAxisLength_LL","MajorAxisLength_LR","MajorAxisLength_UL","MajorAxisLength_UR","Maximum2DDiameterColumn_CR","Maximum2DDiameterColumn_LL","Maximum2DDiameterColumn_LR","Maximum2DDiameterColumn_UL","Maximum2DDiameterColumn_UR","Maximum2DDiameterRow_CR","Maximum2DDiameterRow_LL","Maximum2DDiameterRow_LR","Maximum2DDiameterRow_UL","Maximum2DDiameterRow_UR","Maximum2DDiameterSlice_CR","Maximum2DDiameterSlice_LL","Maximum2DDiameterSlice_LR","Maximum2DDiameterSlice_UL","Maximum2DDiameterSlice_UR","Maximum3DDiameter_CR","Maximum3DDiameter_LL","Maximum3DDiameter_LR","Maximum3DDiameter_UL","Maximum3DDiameter_UR","MeshVolume_CR","MeshVolume_LL","MeshVolume_LR","MeshVolume_UL","MeshVolume_UR","MinorAxisLength_CR","MinorAxisLength_LL","MinorAxisLength_LR","MinorAxisLength_UL","MinorAxisLength_UR","Sphericity_CR","Sphericity_LL","Sphericity_LR","Sphericity_UL","Sphericity_UR","SurfaceArea_CR","SurfaceArea_LL","SurfaceArea_LR","SurfaceArea_UL","SurfaceArea_UR","SurfaceVolumeRatio_CR","SurfaceVolumeRatio_LL","SurfaceVolumeRatio_LR","SurfaceVolumeRatio_UL","SurfaceVolumeRatio_UR","VoxelVolume_CR","VoxelVolume_LL","VoxelVolume_LR","VoxelVolume_UL","VoxelVolume_UR")
# intercept_content$Phenotypes <- names
intercept_content$Intercept <- as.numeric(intercept_content$Intercept)

#输出结果，对应159个表型（少一个minimum_UR因为LDSC遗传度过小无法计算）
writexl::write_xlsx(intercept_content,'LDSC_phenptype_intercept.xlsx')

#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################
#####################################################UR#####################################################



#####################################################宏观#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
rm(list=ls())
install.packages("Hmisc")

library(data.table)
library(Hmisc)

phe=fread("Phenotype.txt")
phe=phe[,-c(1:2)]
sample <- read.csv('35469sample_GWAS.csv')
table(sample$x %in% phe$IID) #T
phe <- subset(phe, IID %in% sample$x)


#提取各个肺叶的值#
UR <- phe[,c(2:33)]


########################UR ########################
#colnames(UR) <- gsub("original_firstorder_|original_shape_|_UR", "", colnames(UR))
colnames(UR) <- gsub("original_|_UR", "", colnames(UR))
UR <- UR[,c(15:32,1:14)] #调整顺序
UR <- UR[,c(3,16,7,12,6,15,17,11,1:2,4:5,8:10,13:14,18,31,19:20,29,21:28,30,32)]


# 计算相关性和p值
cor_results <- rcorr(as.matrix(UR), type = "pearson")
cor_matrix <- cor_results$r
p_matrix <- cor_results$P

# 生成星号矩阵，表示p值的显著性
star_mat <- ifelse(p_matrix < 0.05/32, "**", 
                   ifelse(p_matrix < 0.05, "*", ""))

# 提取宏观矩阵的upper部分为新矩阵的upper
#  cor_matrix 是宏观的原始相关矩阵  
# 创建一个与 cor_matrix 相同大小的新矩阵，并初始化为 NA  
UR_matrix <- matrix(NA, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))  
# 填充上半部分（不包括对角线）  
UR_matrix[upper.tri(UR_matrix, diag = FALSE)] <- cor_matrix[upper.tri(cor_matrix, diag = FALSE)]

# 对星号矩阵进行上半部分的提取
# 创建一个与 cor_matrix 相同大小的新矩阵，并初始化为空字符串
combined_star_mat <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
# 填充上半部分为宏观星号矩阵的上半部分（不包括对角线）
combined_star_mat[upper.tri(combined_star_mat, diag = FALSE)] <- star_mat[upper.tri(star_mat, diag = FALSE)]

#####################################################遗传#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
R
library(ComplexHeatmap)

dir_path <- "./step4_lobe_inner_ldsc/result/UR"

# 获取文件列表
files <- list.files(dir_path, full.names = TRUE)

# 初始化一个空数据框来存储结果
result_df <- data.frame()

# 遍历文件
for (file in files) {
  # 读取文件内容
  content <- readLines(file)
  
  # 判断文件中是否存在 "Genetic Correlation: " 行
  if (!any(grepl("Genetic Correlation: ", content))) {
    next
  }
  
  # 提取 "Genetic Correlation: " 后面到 "(" 之前的内容（不包括冒号后的空格）
  genetic_corr <- grep("Genetic Correlation: ", content, value = TRUE)
  genetic_corr <- sub(".*Genetic Correlation:\\s*(.*?)\\(.*", "\\1", genetic_corr)
  
  # 提取 "P: " 后面到字符串结尾的内容（不包括冒号后的空格）
  p_value <- grep("P: ", content, value = TRUE)
  p_value <- sub(".*P:\\s*", "", p_value)
  
  # 提取文件的 basename 去掉 "_rg.log" 部分作为文件名
  filename <- basename(file)
  filename <- sub("_rg.log$", "", filename)
  
  # 将结果添加到数据框
  result_df <- rbind(result_df, data.frame(File = filename, Genetic_Correlation = genetic_corr, p = p_value))
}

# 在循环结束后拆分文件名
result_df <- transform(result_df, phe1 = sub("\\..*$", "", File), phe2 = sub("^[^.]*\\.", "", File))

# 删除原始文件名列
result_df <- result_df[, c("phe1", "phe2", "Genetic_Correlation", "p")]
result_df$Genetic_Correlation<-as.numeric(result_df$Genetic_Correlation)
result_df$p<-ifelse(result_df$p=="0.",0,result_df$p)
result_df$p<-as.numeric(result_df$p)
result_df$p.adj<-p.adjust(result_df$p,method="fdr")

#########################整理成相关性矩阵####################
# 提取唯一的表型
unique_phe <- unique(c(result_df$phe1, result_df$phe2))

# 创建相关性值矩阵
cor_matrix <- matrix(NA, nrow = length(unique_phe), ncol = length(unique_phe), dimnames = list(unique_phe, unique_phe))
# 去除行名和列名中的 "original_" 和 "_UR"
rownames(cor_matrix) <- sub("original_", "", rownames(cor_matrix))
rownames(cor_matrix) <- sub("_UR", "", rownames(cor_matrix))
colnames(cor_matrix) <- sub("original_", "", colnames(cor_matrix))
colnames(cor_matrix) <- sub("_UR", "", colnames(cor_matrix))

# 创建 p 值矩阵
p_value_matrix <- matrix(NA, nrow = length(unique_phe), ncol = length(unique_phe), dimnames = list(unique_phe, unique_phe))
rownames(p_value_matrix) <- sub("original_", "", rownames(p_value_matrix))
rownames(p_value_matrix) <- sub("_UR", "", rownames(p_value_matrix))
colnames(p_value_matrix) <- sub("original_", "", colnames(p_value_matrix))
colnames(p_value_matrix) <- sub("_UR", "", colnames(p_value_matrix))

# 创建 p.adj 值矩阵
p_adj_matrix <- matrix(NA, nrow = length(unique_phe), ncol = length(unique_phe), dimnames = list(unique_phe, unique_phe))
rownames(p_adj_matrix) <- sub("original_", "", rownames(p_adj_matrix))
rownames(p_adj_matrix) <- sub("_UR", "", rownames(p_adj_matrix))
colnames(p_adj_matrix) <- sub("original_", "", colnames(p_adj_matrix))
colnames(p_adj_matrix) <- sub("_UR", "", colnames(p_adj_matrix))

# 填充相关性值矩阵和 p 值矩阵
for (i in 1:nrow(result_df)) {
  row_index <- which(unique_phe == result_df[i, "phe1"])
  col_index <- which(unique_phe == result_df[i, "phe2"])
  
  cor_matrix[row_index, col_index] <- result_df[i, "Genetic_Correlation"]
  cor_matrix[col_index, row_index] <- result_df[i, "Genetic_Correlation"]
  
  p_value_matrix[row_index, col_index] <- result_df[i, "p"]
  p_value_matrix[col_index, row_index] <- result_df[i, "p"]

  p_adj_matrix[row_index, col_index] <- result_df[i, "p.adj"]
  p_adj_matrix[col_index, row_index] <- result_df[i, "p.adj"]
}

# 对角线元素设置为 1
diag(cor_matrix) <- 1
diag(p_value_matrix) <- 1
diag(p_adj_matrix) <- 1

cor_matrix[cor_matrix > 1] <- 1
cor_matrix[cor_matrix < -1] <- -1

# 定义目标顺序
target_order <- c("firstorder_Energy", "firstorder_TotalEnergy","firstorder_Maximum", "firstorder_Range",
                  "firstorder_Kurtosis", "firstorder_Skewness",
                  "firstorder_Uniformity", "firstorder_10Percentile",
                  "firstorder_90Percentile", "firstorder_Entropy", "firstorder_InterquartileRange",
                  "firstorder_MeanAbsoluteDeviation", "firstorder_Mean", "firstorder_Median",
                  "firstorder_RobustMeanAbsoluteDeviation", "firstorder_RootMeanSquared",
                  "firstorder_Variance", "shape_SurfaceVolumeRatio", "shape_Elongation",
                  "shape_Flatness", "shape_Sphericity", "shape_LeastAxisLength",
                  "shape_MajorAxisLength", "shape_Maximum2DDiameterColumn",
                  "shape_Maximum2DDiameterRow", "shape_Maximum2DDiameterSlice",
                  "shape_Maximum3DDiameter", "shape_MeshVolume", "shape_MinorAxisLength",
                  "shape_SurfaceArea", "shape_VoxelVolume")     

# 按照目标顺序重新排列矩阵
cor_matrix <- cor_matrix[target_order, target_order]
p_value_matrix <- p_value_matrix[target_order, target_order]
p_adj_matrix <- p_adj_matrix[target_order, target_order]

# 创建一个矩阵来存储带有星号的矩阵
star_mat <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
# 找出 p 值小于 0.05 的位置，并在 star_mat 中相应位置填入星号
star_mat[p_value_matrix < 0.05/31] <- "*"
# 将 cor_matrix 的行名和列名赋给 star_mat
rownames(star_mat) <- rownames(cor_matrix)
colnames(star_mat) <- colnames(cor_matrix)


### ### ### ### ### 遗传矩阵缺失了minimum，需要将其补进去
###数据
ur_names <- rownames(UR_matrix)
cor_names <- rownames(cor_matrix)

# 找出 cor_matrix 中缺失的变量
missing_variable <- setdiff(ur_names, cor_names)
missing_variable
#[1] "firstorder_Minimum"
# 创建一个新的 32x32 矩阵，初始为 NA
new_cor_matrix <- matrix(NA, nrow = 32, ncol = 32)
# 将 cor_matrix 的原有数据放在新矩阵的对应位置
new_cor_matrix[1:31, 1:31] <- cor_matrix
# 设置新矩阵的行名和列名
rownames(new_cor_matrix) <- c(cor_names, missing_variable)
colnames(new_cor_matrix) <- c(cor_names, missing_variable)
# 按照 UR_matrix 的顺序重新排列行和列
new_cor_matrix <- new_cor_matrix[ur_names, ur_names]
###星号
ur_names <- rownames(combined_star_mat)
star_names <- rownames(star_mat)
# 找出 star_mat 中缺失的变量
missing_variable <- setdiff(ur_names, star_names)
missing_variable
#[1] "firstorder_Minimum"
# 创建一个新的 32x32 矩阵，初始为 空字符串
new_star_mat <- matrix("", nrow = 32, ncol = 32)
# 将 star_mat 的原有数据放在新矩阵的对应位置
new_star_mat[1:31, 1:31] <- star_mat
# 设置新矩阵的行名和列名
rownames(new_star_mat) <- c(star_names, missing_variable)
colnames(new_star_mat) <- c(star_names, missing_variable)
# 按照 combined_star_mat 的顺序重新排列行和列
new_star_mat <- new_star_mat[ur_names, ur_names]

# 提取遗传矩阵的lower部分为新矩阵的lower
# new_cor_matrix 是遗传的原始相关矩阵  
# 填充下半部分（不包括对角线）  
UR_matrix[lower.tri(UR_matrix, diag = FALSE)] <- new_cor_matrix[lower.tri(new_cor_matrix, diag = FALSE)]

# 对星号矩阵进行下半部分的提取
# 填充下半部分为遗传星号矩阵的下半部分（不包括对角线）
combined_star_mat[lower.tri(combined_star_mat, diag = FALSE)] <- new_star_mat[lower.tri(new_star_mat, diag = FALSE)]

### ### ### ### ###画图
pdf("Heatmap_UR.pdf", width = 8, height = 8)

pheatmap(UR_matrix,name="rg_UR",
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 8, cellheight = 8,  #设置单元格大小
         scale="none",
         display_numbers = combined_star_mat, #设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), #设置颜色 #我设置的好丑
         annotation_legend="rg",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("90"),  #设置列标签的倾斜角度
         show_rownames=TRUE
)

dev.off()

#############same for the other 4 lobes#############
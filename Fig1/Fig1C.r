R
library(ComplexHeatmap)

load("UR_matrix_all.RData") #UR_matrix
load("UL_matrix_all.RData") #UL_matrix
load("CR_matrix_all.RData") #CR_matrix
load("LR_matrix_all.RData") #LR_matrix
load("LL_matrix_all.RData") #LL_matrix

matrices <- list(UR_matrix, UL_matrix, CR_matrix, LR_matrix, LL_matrix)

# 获取矩阵的维度（假设所有矩阵的维度相同）
nrow_mat <- nrow(UR_matrix)
ncol_mat <- ncol(UR_matrix)

# 初始化一个矩阵用于存储结果
mean_matrix <- matrix(NA, nrow = nrow_mat, ncol = ncol_mat)

# 遍历每个格子计算均值
for (i in 1:nrow_mat) {
  for (j in 1:ncol_mat) {
    # 提取当前格子在5个矩阵中的值
    values <- sapply(matrices, function(mat) mat[i, j])
    mean_matrix[i, j] <- mean(values, na.rm = TRUE)
  }
}

# 将UR_matrix的行名和列名赋给mean_matrix
rownames(mean_matrix) <- rownames(UR_matrix)
colnames(mean_matrix) <- colnames(UR_matrix)

### ### ### ### ###plot
pdf("Heatmap_all.pdf", width = 8, height = 8)

pheatmap(mean_matrix,name="rg_all",
         cluster_rows=FALSE, 
         cluster_cols = FALSE, 
         border = "white" ,
         cellwidth = 8, cellheight = 8,  
         scale="none",
         fontsize_number = 6,
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), 
         annotation_legend="rg",
         fontsize_row = 6, 
         fontsize_col = 6,
         angle_col = c("90"),
         show_rownames=TRUE
)

dev.off()
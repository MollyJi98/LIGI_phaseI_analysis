###################TF,GF与表型的关系绘图###################
library(openxlsx)
library(ComplexHeatmap)
library(extrafont)

data <- read.xlsx('TF_Results_ALL_for_plot.xlsx')
firstorder <- read.xlsx('TF_Results_ALL_for_plot.xlsx',sheet = 2)
shape <- read.xlsx('TF_Results_ALL_for_plot.xlsx',sheet = 3)


###################相关性热图-与first order表型###################
dim(firstorder)
data_h <- firstorder[,c(1,7:96)]
row.names(data_h) <- data_h$Marker
data_h <- data_h[,-1]

# 提取列名
columns <- colnames(data_h)
# 自定义排序顺序
order_suffix <- c("UR", "CR", "LR", "UL", "LL")
# 提取尾缀
suffixes <- sapply(columns, function(x) sub(".*_", "", x))
# 根据自定义顺序排序列名
sorted_columns <- columns[order(match(suffixes, order_suffix))]
# 重新排列数据框的列
data_h <- data_h[, sorted_columns]

#去除列尾的标签
new_colnames <- sapply(sorted_columns, function(x) sub("_.*", "", x))
colnames(data_h) <- new_colnames

data_heatmap <- as.matrix(data_h) #6,90

# 创建一个矩阵来存储带有星号的矩阵
star_mat <- matrix("", nrow = nrow(data_heatmap), ncol = ncol(data_heatmap))
# 找出 p 值小于 0.05 的位置，并在 star_mat 中相应位置填入星号
star_mat[data_heatmap < 0.05] <- "*"
star_mat[data_heatmap < 0.05/90] <- "**"
#star_mat <- t(star_mat)

pdf("Heatmap_TF_GF_firstorder_assoc-v2.pdf", width = 15, height = 10)

pheatmap(-log(data_heatmap),
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 7, cellheight = 7,  #设置单元格大小
         scale="none",
         display_numbers = star_mat, #设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "white","firebrick3" ))(50), #设置颜色
         annotation_legend="rg",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("45"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)

dev.off()

###################相关性热图-与shape表型###################
dim(shape)
data_h <- shape[,c(1,7:76)]
row.names(data_h) <- data_h$Marker
data_h <- data_h[,-1]

# 提取列名
columns <- colnames(data_h)
# 自定义排序顺序
order_suffix <- c("UR", "CR", "LR", "UL", "LL")
# 提取尾缀
suffixes <- sapply(columns, function(x) sub(".*_", "", x))
# 根据自定义顺序排序列名
sorted_columns <- columns[order(match(suffixes, order_suffix))]
# 重新排列数据框的列
data_h <- data_h[, sorted_columns]

#去除列尾的标签
new_colnames <- sapply(sorted_columns, function(x) sub("_.*", "", x))
colnames(data_h) <- new_colnames

data_heatmap <- as.matrix(data_h) #33,70

# 创建一个矩阵来存储带有星号的矩阵
star_mat <- matrix("", nrow = nrow(data_heatmap), ncol = ncol(data_heatmap))
# 找出 p 值小于 0.05 的位置，并在 star_mat 中相应位置填入星号
star_mat[data_heatmap < 0.05] <- "*"
star_mat[data_heatmap < 0.05/70] <- "**"
#star_mat <- t(star_mat)

pdf("Heatmap_TF_GF_shape_assoc-v2.pdf", width = 15, height = 10)

pheatmap(-log(data_heatmap),
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 7, cellheight = 7,  #设置单元格大小
         scale="none",
         display_numbers = star_mat, #设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "white","firebrick3" ))(50), #设置颜色
         annotation_legend="rg",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("45"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)

dev.off()

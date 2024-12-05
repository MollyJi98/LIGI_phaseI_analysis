###############################################################################################################
########################################绘制ldsc的热图(指定顺序版)################################################
###############################################################################################################

library(openxlsx)
library(ComplexHeatmap)


ldsc <- read.xlsx("./LDSC_res_for_plot.xlsx",sheet = 1,colNames = T,rowNames = T)
ldsc_P <- read.xlsx("./LDSC_res_for_plot.xlsx",sheet = 2,colNames = T,rowNames = T)
ldsc_P.adj <- read.xlsx("./LDSC_res_for_plot.xlsx",sheet = 3,colNames = T,rowNames = T)

#first order#
first <- ldsc[c(1:10,16:25,36:55,61:70,86:90,26:30,71:75,81:85,11:15,76:80,31:35,56:60),c(1,2,3,4,5,6,7)] 
first_P<- ldsc_P[c(1:10,16:25,36:55,61:70,86:90,26:30,71:75,81:85,11:15,76:80,31:35,56:60),c(1,2,3,4,5,6,7)]
first_P.adj<- ldsc_P.adj[c(1:10,16:25,36:55,61:70,86:90,26:30,71:75,81:85,11:15,76:80,31:35,56:60),c(1,2,3,4,5,6,7)]


# 创建一个矩阵来存储带有星号的矩阵
star_mat <- matrix("", nrow = nrow(first), ncol = ncol(first))
# 找出 p 值小于 0.05 的位置，并在 star_mat 中相应位置填入星号
star_mat[first_P < 0.05] <- "*"
star_mat[first_P.adj < 0.05] <- "**"
star_mat <- t(star_mat)

pdf("Heatmap_LDSC_firstorder_rerange.pdf", width = 12, height = 6)

pheatmap(as.matrix(t(first)),
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 8, cellheight = 8,  #设置单元格大小
         scale="none",
         display_numbers = star_mat, #设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), #设置颜色 #我设置的好丑
         annotation_legend="rg",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("90"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)

dev.off()


#shape#
shape <- ldsc[c(101:140,146:150,156:160,91:100,141:145,151:155),c(1,2,3,4,5,6,7)]
shape_P<- ldsc_P[c(101:140,146:150,156:160,91:100,141:145,151:155),c(1,2,3,4,5,6,7)]
shape_P.adj<- ldsc_P.adj[c(101:140,146:150,156:160,91:100,141:145,151:155),c(1,2,3,4,5,6,7)]


# 创建一个矩阵来存储带有星号的矩阵
star_mat <- matrix("", nrow = nrow(shape), ncol = ncol(shape))
# 找出 p 值小于 0.05 的位置，并在 star_mat 中相应位置填入星号
star_mat[shape_P < 0.05] <- "*"
star_mat[shape_P.adj < 0.05] <- "**"
star_mat <- t(star_mat)

pdf("Heatmap_LDSC_shape_rerange.pdf", width = 12, height = 6)

pheatmap(as.matrix(t(shape)),
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 8, cellheight = 8,  #设置单元格大小
         scale="none",
         display_numbers = star_mat, #设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), #设置颜色 #我设置的好丑
         annotation_legend="rg",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("90"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)

dev.off()

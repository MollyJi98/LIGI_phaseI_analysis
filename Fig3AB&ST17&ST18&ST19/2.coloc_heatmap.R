###################################################################################################################
########################################共定位结果画图(不考虑肺叶，画两类总体的图)########################################
###################################################################################################################

rm(list=ls())
library(ComplexHeatmap)
library(reshape2)
library(gtools)
library(tidyr)  
 
coloc <- read.csv('coloc_results_for_use.csv') #挑选的显著的pph4>0.7且表型对应区域达到显著的共定位-表型对，剩下的pph4>0.7,且没有报道的共定位基因！！！！！！都是ST19的内容
all=read.csv('allegene_coloc_result.csv') #all coloc result
head(all)


########################################提取pph4结果（firstorder）########################################
#提取17个不重复的firstorder的共定位基因
firstorder=read.xlsx("coloc_gene_for_plot-jc.xlsx",sheet=1) #firstorder所有要画图的基因

# 提取所有结果中对应基因的行
filtered <- all[grep("^original_firstorder", all$phenotype), ]

# 保留gene_id在firstorder$gene_id中的行，保留画图所需要的三列变量
filtered <- filtered[filtered$gene_id %in% firstorder$gene_id, ]
filtered$reg_gene <- paste(filtered$Cytoband,'(',filtered$gene_name,')',sep = "")
filtered <- filtered[,c('phenotype','PPH4','reg_gene')]

# 使用 dcast 函数将 filtered 数据框变换为适合绘制热图的矩阵
heatmap_matrix <- dcast(filtered, phenotype ~ reg_gene, value.var = "PPH4")
heatmap_matrix2 <- heatmap_matrix %>%  
  separate(phenotype, into = c("part1", "part2", "part3", "part4"), sep = "_") %>%  
  mutate(phenotype_new = paste(part1, part2, part3, sep = "_")) %>%  
  select(-part1, -part2, -part3, -part4)  # 删除拆分后的列（如果需要）  

heatmap_matrix3 <- heatmap_matrix2 %>%  
  group_by(phenotype_new) %>%  
  summarise(across(1:17, ~max(., na.rm = TRUE)), .groups = 'drop')   %>% as.data.frame()

# 将 phenotype 列设置为行名
rownames(heatmap_matrix3) <- heatmap_matrix3$phenotype_new
heatmap_matrix3$phenotype_new <- NULL  # 移除 phenotype 列

# 将矩阵转置，以便 gene_name 列作为列名
heatmap_matrix3 <- t(heatmap_matrix3) #17  90
# 获取列名
colnames(heatmap_matrix3) <- gsub("original_firstorder_", "", colnames(heatmap_matrix3))

# 提取行名
rownames <- rownames(heatmap_matrix3)
# 定义一个函数，根据行名的两个部分进行排序
custom_sort <- function(rownames) {
  numeric_part <- as.numeric(gsub("[^0-9].*$", "", rownames))
  non_numeric_part <- gsub("^[0-9]*", "", rownames)
  sorted_indices <- order(numeric_part, non_numeric_part)
  return(sorted_indices)
}

sorted_indices <- custom_sort(rownames)
sorted_heatmap_matrix <- heatmap_matrix3[sorted_indices, ] #17,19


# 指定要保存的 PDF 文件名
pdf("Heatmap_firstorder_colocsig_0.7_combinelobe.pdf", width = 12, height = 6)

pheatmap(sorted_heatmap_matrix,
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 8, cellheight = 8,  #设置单元格大小
         display_numbers = matrix(ifelse(sorted_heatmap_matrix >= 0.9, "**",ifelse(sorted_heatmap_matrix >= 0.7, "*","")), nrow = nrow(sorted_heatmap_matrix)),#设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), #设置颜色 #我设置的好丑
         annotation_legend="pph4",
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("45"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)


dev.off()



########################################提取pph4结果（shape）########################################
#提取28个不重复的firstorder的共定位基因
shape=read.xlsx("coloc_gene_for_plot-jc.xlsx",sheet=2) #

# 提取所有结果中对应基因的行
filtered <- all[grep("^original_shape", all$phenotype), ]

# 保留gene_id在firstorder$gene_id中的行，保留画图所需要的三列变量
filtered <- filtered[filtered$gene_id %in% shape$gene_id, ]
filtered$reg_gene <- paste(filtered$Cytoband,'(',filtered$gene_name,')',sep = "")
filtered <- filtered[,c('phenotype','PPH4','reg_gene')]

# 使用 dcast 函数将 filtered 数据框变换为适合绘制热图的矩阵
heatmap_matrix <- dcast(filtered, phenotype ~ reg_gene, value.var = "PPH4")
heatmap_matrix2 <- heatmap_matrix %>%  
  separate(phenotype, into = c("part1", "part2", "part3", "part4"), sep = "_") %>%  
  mutate(phenotype_new = paste(part1, part2, part3, sep = "_")) %>%  
  select(-part1, -part2, -part3, -part4)  # 删除拆分后的列（如果需要）  

heatmap_matrix3 <- heatmap_matrix2 %>%  
  group_by(phenotype_new) %>%  
  summarise(across(1:28, ~max(., na.rm = TRUE)), .groups = 'drop')   %>% as.data.frame()

# 将 phenotype 列设置为行名
rownames(heatmap_matrix3) <- heatmap_matrix3$phenotype_new
heatmap_matrix3$phenotype_new <- NULL  # 移除 phenotype 列

# 将矩阵转置，以便 gene_name 列作为列名
heatmap_matrix3 <- t(heatmap_matrix3) #28  14
# 获取列名
colnames(heatmap_matrix3) <- gsub("original_shape_", "", colnames(heatmap_matrix3))


# 提取行名
rownames <- rownames(heatmap_matrix3)
# 定义一个函数，根据行名的两个部分进行排序
custom_sort <- function(rownames) {
  numeric_part <- as.numeric(gsub("[^0-9].*$", "", rownames))
  non_numeric_part <- gsub("^[0-9]*", "", rownames)
  sorted_indices <- order(numeric_part, non_numeric_part)
  return(sorted_indices)
}

sorted_indices <- custom_sort(rownames)
sorted_heatmap_matrix <- heatmap_matrix3[sorted_indices, ] 


pdf("Heatmap_shape_colocsig_0.7_combinelobe.pdf", width = 12, height = 6)

pheatmap(sorted_heatmap_matrix,
         cluster_rows=FALSE, #是否对行进行聚类
         cluster_cols = FALSE, #是否对列进行聚类
         border = "white" ,#设置边框颜色，删除边框直接F
         cellwidth = 8, cellheight = 8,  #设置单元格大小
         display_numbers = matrix(ifelse(sorted_heatmap_matrix >= 0.9, "**",ifelse(sorted_heatmap_matrix >= 0.7, "*","")), nrow = nrow(sorted_heatmap_matrix)),#设置单元格数值
         # number_color = "#4daf4a",# 设置颜色
         fontsize_number = 6, # 设置数值字体大小   
         color = colorRampPalette(c( "#5b98e0","white","firebrick3" ))(50), #设置颜色 #我设置的好丑
         fontsize_row = 6, #设置标签的大小
         fontsize_col = 6,
         angle_col = c("45"),  #设置列标签的倾斜角度
         show_rownames=TRUE,
         show_colnames=TRUE
)

dev.off()


##################################完成后在ai里拼图##################################
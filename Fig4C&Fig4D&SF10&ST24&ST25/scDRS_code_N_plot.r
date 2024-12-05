#########################h5ad文件和cov文件准备#########################
library(Seurat)
library(SeuratDisk)

load('Cell_FL_scDRS.Rdata') #单细胞数据
ls() #查看数据集

table(All.merge$stage)
#    5  6.86     9    11    15    18    20    22 
# 4203  4790  8280  9531 12071 13170  7490 12217 

da_15_22 <- subset(All.merge, subset = stage >= 15)

save(da_15_22,file="Cell_FL_scDRS_15_22.Rdata") #26568 44948

SaveH5Seurat(da_15_22, filename = "Cell_FL_scDRS_15_22.h5Seurat",overwrite = T)
Convert("Cell_FL_scDRS_15_22.h5Seurat", dest = "h5ad",assay="RNA",overwrite = T)

# 协变量文件提取
covariate_data <- da_15_22@meta.data

# 假设我们只需要 batch 和 gender 作为协变量
covariate_data_sorted <- as.data.frame(covariate_data[, c("gender","nCount_RNA","nFeature_RNA")])
rownames(covariate_data_sorted)<-rownames(covariate_data)
names(covariate_data_sorted)[1]<-"gender"
covariate_data_sorted$gender<-ifelse(covariate_data_sorted$gender=="F",2,1)
covariate_data_sorted$Cell_ID<-rownames(covariate_data_sorted) #注意一定要重新定义一列！！！！和h5ad文件名一致！！！！！
covariate_data_sorted<-covariate_data_sorted[,c("Cell_ID","gender","nCount_RNA","nFeature_RNA")]

# 检查协变量数据的细胞ID与h5ad文件中的细胞ID是否一致
cov_cells <- rownames(covariate_data_sorted)
h5ad_cells <- colnames(da_15_22)  # h5ad文件中的细胞名称

# 保存为 scdrs 协变量文件
write.table(covariate_data_sorted, file = "covariates_for_scdrs.cov", row.names = FALSE,quote=F,sep="\t")

# 绘制UMAP图(broad_celltype)
# 将缺失值替换为 "NK"
da_15_22@meta.data$broad_celltype[is.na(da_15_22@meta.data$broad_celltype)] <- "NK"

# 将 broad_celltype 转换为 factor 类型，并设置所需的 levels
da_15_22@meta.data$broad_celltype <- factor(
  da_15_22@meta.data$broad_celltype,
  levels = c("B", "Chondrocyte", "Distal epithelial", "Fibroblast", "Lymph endothelial", 
             "Meg-ery", "Mesothelial", "Myofibro & SMC", "Other myeloid", "PNS", 
             "Proximal epithelial", "T & ILC", "Vas endothelial", "NK"))

da_15_22@meta.data$new_celltype<-factor(da_15_22@meta.data$new_celltype)
unique_celltypes <- length(unique(da_15_22@meta.data$new_celltype)) #查看new_celltype的种类数
print(unique_celltypes)


mycolor <- rainbow(14, s = 0.6, v = 0.8)
pdf(file='broad_celltype_anno.pdf',width=10,heigh=6)
DimPlot(object=da_15_22,pt.size=0.01,group.by='broad_celltype',raster=FALSE,cols=mycolor, label = T,repel = T)
dev.off()



##########################################################################################################
########################################all magma结果，同一个表型合并########################################
##########################################################################################################

scdrs compute_score \
    --h5ad-file Cell_FL_scDRS_15_22.h5ad \
    --h5ad-species human \
    --gs-file scdrs_gene_set_file_all_magma_top1000gene_combined.tsv \
    --gs-species human \
    --cov-file covariates_for_scdrs.cov \
    --flag-filter-data True \
    --flag-raw-count False \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder /top_1000_magma/result/




# 批量处理scdrs结果
library(data.table)
library(dplyr)
library(writexl)

# 读取covariate_data文件
load('Cell_FL_scDRS.Rdata')
da_15_22 <- subset(All.merge, subset = stage >= 15)
covariate_data <- da_15_22@meta.data
covariate_data$V1 <- rownames(covariate_data)

# 获取所有以.score.gz结尾的文件
files <- list.files("./top_1000_magma/result", pattern = "\\.score\\.gz$", full.names = TRUE)

# 初始化一个空列表用于存储每个文件的结果
all_results <- list()

# 遍历每个文件
for (file in files) {
    # 获取文件名，去掉路径和扩展名，作为GWAS_pheno
    GWAS_pheno <- gsub(".score.gz", "", basename(file))
    # 读取数据并与covariate_data合并
    first_res <- fread(file)
    first_res_combind <- merge(first_res, covariate_data, by = "V1")
    # 对pval进行fdr校正
    first_res_combind$pval_adj <- p.adjust(first_res_combind$pval, method = "fdr")
    # 计算 broad_celltype 的各类比例
    result <- first_res_combind %>%
      group_by(broad_celltype) %>%
      summarize(
        total_count = n(),  # 每种 broad_celltype 的总数
        count_pval_lt_0_05 = sum(pval < 0.05),  # pval < 0.05 的数量
        count_pval_lt_0_1_adj = sum(pval_adj < 0.1),  # pval_adj < 0.1 的数量
        count_pval_lt_0_05_adj = sum(pval_adj < 0.05),  # pval_adj < 0.05 的数量
        proportion_pval_lt_0_05 = count_pval_lt_0_05 / total_count,  # pval < 0.05 的比例
        proportion_pval_lt_0_1_adj = count_pval_lt_0_1_adj / total_count,  # pval_adj < 0.1 的比例
        proportion_pval_lt_0_05_adj = count_pval_lt_0_05_adj / total_count  # pval_adj < 0.05 的比例
      ) %>% 
      mutate(GWAS_pheno = GWAS_pheno) %>%  # 增加GWAS_pheno列
      as.data.frame()
    # 将结果添加到列表中
    all_results[[GWAS_pheno]] <- result
}

# 将所有结果合并成一个大数据框
final_result <- bind_rows(all_results)

# 输出结果到Excel
write.table(final_result, "combine_res_all_files_propotion.txt",quote=F,sep="\t",row.name=F)





# scdrs-downstream
input_dir="./top_1000_magma/result"
h5ad_file="Cell_FL_scDRS_15_22.h5ad"
output_dir="./top_1000_magma/downstream"

# 循环遍历所有 .full_score.gz 文件
for score_file in "$input_dir"/*.full_score.gz; do
    # 获取文件的名称（去除路径和扩展名）
    base_name=$(basename "$score_file" .full_score.gz)
    
    # 运行 scDRS 下游处理
    scdrs perform-downstream \
    --h5ad-file "$h5ad_file" \
    --score-file "$score_file" \
    --out-folder "$output_dir" \
    --group-analysis broad_celltype \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count True
done

# 合并downstream和比例的结果
library(data.table)
library(openxlsx)
library(ggplot2)
library(dplyr)

propotion<-fread('combine_res_all_files_propotion.txt')

# 合并所有表型mctest的结果
input_dir <- "./top_1000_magma/downstream"
output_file <- "combined_result_mctest_top_1000_magma.txt"

# 获取所有文件名
file_list <- list.files(path = input_dir, pattern = "\\.scdrs_group\\.broad_celltype$", full.names = TRUE)
# 初始化一个空的数据框用于合并
combined_data <- data.frame()

# 遍历每个文件
for (file in file_list) {
  # 读取文件内容
  temp_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # 提取文件名去掉路径和扩展名
  file_name <- basename(file)
  gwas_pheno <- sub("\\.scdrs_group\\.broad_celltype$", "", file_name)
  # 添加新的一列 GWAS_pheno
  temp_data$GWAS_pheno <- gwas_pheno
  temp_data$assoc_mcp_adj<-p.adjust(temp_data$assoc_mcp,method="fdr")
  temp_data$hetero_mcp_adj<-p.adjust(temp_data$hetero_mcp,method="fdr")
  # 合并到总的数据框
  combined_data <- rbind(combined_data, temp_data)
}


mc<-combined_data # 也就是计算完的mc结果，接下来和propotion合并
final_result <- data.frame()

# 获取 propotion 和 mc 数据框中共有的 GWAS_pheno 值
gwas_pheno_list <- intersect(unique(propotion$GWAS_pheno), unique(mc$GWAS_pheno))

# 遍历每个 GWAS_pheno
for (pheno in gwas_pheno_list) {
  # 筛选当前 GWAS_pheno 的数据
  propotion_subset <- propotion[propotion$GWAS_pheno == pheno, ]
  mc_subset <- mc[mc$GWAS_pheno == pheno, ]
  # 合并两个子集
  merged_subset <- merge(propotion_subset[,c(1,3,6)], mc_subset, by.x = "broad_celltype", by.y = "group")
  # 将合并结果追加到 final_result
  final_result <- rbind(final_result, merged_subset)
}

# 如果需要保存结果到文件
write.table(final_result, "propotion_mctest_all_phe_res.txt", row.names = FALSE,sep="\t",quote=F)


##################基于结果绘制热图heatmap_plot####################
library(ggplot2)
library(reshape2)
library(dplyr)


# 数据预处理
heatmap_data <- final_result %>%
  select(GWAS_pheno, broad_celltype, proportion_pval_lt_0_05, assoc_mcp_adj, assoc_mcp) %>%
  mutate(
    GWAS_pheno = gsub("original_firstorder_|original_shape_", "", GWAS_pheno),  # 去除前缀
    proportion_bin = cut(
      proportion_pval_lt_0_05,
      breaks = seq(0, 1, by = 0.1),
      labels = paste0(seq(0, 90, by = 10), "%-", seq(10, 100, by = 10), "%"),
      include.lowest = TRUE
    ),
    significance = case_when(
      assoc_mcp_adj < 0.05 ~ "**",  
      assoc_mcp < 0.05 ~ "*",      
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(GWAS_pheno = factor(GWAS_pheno, levels = rev(unique(GWAS_pheno)))) 

# 定义颜色
color_palette <- colorRampPalette(c("#fff9f9", "firebrick3"))(10)

# 定义细胞类型顺序
cell_order <- c(
  # 白细胞和红细胞(leukocyte/erythroid)
  "B", "NK", "T & ILC", "Other myeloid", "Meg-ery",
  # epithelial
  "Distal epithelial", "Proximal epithelial",
  # endothelial
  "Lymph endothelial", "Vas endothelial",
  # fibroblast
  "Fibroblast", "Myofibro & SMC",
  # other
  "Mesothelial", "Chondrocyte", "PNS"
)
heatmap_data$broad_celltype <- factor(heatmap_data$broad_celltype, levels = cell_order)

# 分割数据
firstorder_data <- heatmap_data %>% slice(1:252)  # 前252行(18*14)
shape_data <- heatmap_data %>% slice(253:448)     # 后196行(14*14)

# 绘图函数
plot_heatmap <- function(data, title, output_file) {
  p <- ggplot(data, aes(x = broad_celltype, y = GWAS_pheno, fill = proportion_bin)) +
    geom_tile(color = "white") +  # 绘制热图方块
    scale_fill_manual(
      values = color_palette,
      na.value = "gray90",
      guide = guide_legend(title = "Proportion")
    ) +
    theme_minimal() +
    labs(x = "Broad Cell Type", y = "GWAS Phenotype", title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),  # X轴字体黑色
      axis.text.y = element_text(size = 6, color = "black"),  # Y轴字体黑色
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),  # 标题居中且加大字号
      axis.ticks = element_line(color = "black", size = 0.5),  # 刻度线黑色
      panel.grid.major = element_blank(),  # 移除主网格线
      panel.grid.minor = element_blank()  # 移除次网格线
    ) +
    coord_fixed(ratio = 0.9) +
    # 添加显著性标记
    geom_text(data = data %>% filter(!is.na(significance)),
              aes(label = significance),
              color = "black", size = 3, vjust = 0.5, hjust = 0.5)

  # 保存热图为 PDF 文件
  ggsave(output_file, plot = p, width = 6, height = 5)
}

# 输出图表
plot_heatmap(firstorder_data, "Firstorder Phenotypes", "firstorder_heatmap.pdf")
plot_heatmap(shape_data, "Shape Phenotypes", "shape_heatmap.pdf")

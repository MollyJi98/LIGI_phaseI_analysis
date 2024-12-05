##########对于10traits_gwas_catalog_reported_merge_effect_noNA的位点，按照Phenotype-traitname-rsid进行去重后，再次处理##########
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

catalog <- read.xlsx('10traits_gwas_catalog_reported_merge_effect_noNA.xlsx')

# 定义一个函数来计算每一行从第七列开始的奇数列 < 0.05 的个数（对应的都是p值列）
count_odd_columns_lt_0.05 <- function(row) {
  odd_columns <- as.numeric(row[seq(7, length(row), by = 2)]) # 确保转换为数值类型
  sum(odd_columns < 0.05, na.rm = TRUE)
}

# 定义一个函数来计算每一行从第七列开始的奇数列 < 1e-4 的个数
count_odd_columns_lt_1e_4 <- function(row) {
  odd_columns <- as.numeric(row[seq(7, length(row), by = 2)]) # 确保转换为数值类型
  sum(odd_columns < 1e-4, na.rm = TRUE)
}

# 定义一个函数来计算每一行从第七列开始的奇数列 < 5e-8 的个数
count_odd_columns_lt_5e_8 <- function(row) {
  odd_columns <- as.numeric(row[seq(7, length(row), by = 2)]) # 确保转换为数值类型
  sum(odd_columns < 5e-8, na.rm = TRUE)
}

# 定义一个函数来计算每一行从第七列开始的奇数列 < 0.05/160 的个数
count_odd_columns_lt_adj <- function(row) {
  odd_columns <- as.numeric(row[seq(7, length(row), by = 2)]) # 确保转换为数值类型
  sum(odd_columns < 0.05/160, na.rm = TRUE)
}

# 对于每一行计算 < 0.05 的个数，并添加到新列 count_lt_0.05
catalog <- catalog %>%
  mutate(count_lt_0.05 = apply(., 1, count_odd_columns_lt_0.05),
         count_lt_1e_4 = apply(., 1, count_odd_columns_lt_1e_4),
         count_lt_5e_8 = apply(., 1, count_odd_columns_lt_5e_8),
         count_lt_adj = apply(., 1, count_odd_columns_lt_adj))


# 统计不同 Phenotype 中 count_lt_0.05 > 0 的个数和总行数，并计算比例
result <- catalog %>%
  group_by(Phenotype) %>%
  summarise(
    count_gt_0_5 = sum(count_lt_0.05 > 0, na.rm = TRUE),  # 计算 count_lt_0.05 > 0 的个数
    count_gt_1e_4 = sum(count_lt_1e_4 > 0, na.rm = TRUE),
    count_gt_5e_8 = sum(count_lt_5e_8 > 0, na.rm = TRUE),
    count_lt_adj = sum(count_lt_adj > 0, na.rm = TRUE),
    total_rows = n(),  # 计算每个 Phenotype 的总行数
    proportion_0_5 = count_gt_0_5 / total_rows,  # 计算比例
    proportion_1e_4 = count_gt_1e_4 / total_rows,
    proportion_5e_8 = count_gt_5e_8 / total_rows,
    proportion_adj = count_lt_adj / total_rows
  )

write.table(catalog,'10traits_gwas_catalog_reported_merge_effect_noNA_count.xlsx',quote = F,sep = '\t',row.names = F)
write.table(result,'10traits_gwas_catalog_reported_merge_effect_noNA_count_porpotion_0.05.xlsx',quote = F,sep = '\t',row.names = F)


############################画比例图############################
result <- as.data.frame(result)

# 按照proportion排序
result <- result %>%
  arrange(total_rows)

result <- result[order(result$total_rows, decreasing = TRUE), ] #调整顺序，从高到低
result$Phenotype <- factor(result$Phenotype, levels = result$Phenotype)

# 绘制柱状图
p1 <- ggplot(result[c(1:7,10),], aes(x = Phenotype)) +
  # 绘制 total_rows 的柱状图，设定 alpha 值使其半透明
  geom_bar(aes(y = total_rows), stat = "identity", fill = "lightgray", alpha = 0.6, width = 0.8) + 
  # 绘制 count_gt_0_5 的柱状图，设定不同的颜色并保持半透明
  geom_bar(aes(y = count_gt_0_5), stat = "identity", fill = "#5b98e0", alpha = 0.9, width = 0.8) +
  # 绘制 count_lt_adj 的柱状图，设定不同的颜色并保持半透明
  geom_bar(aes(y = count_lt_adj), stat = "identity", fill = "#ffa500", alpha = 0.9, width = 0.8) +
  # 绘制 count_gt_5e_8 的柱状图，设定不同的颜色并保持半透明
  geom_bar(aes(y = count_gt_5e_8), stat = "identity", fill = "firebrick3", alpha = 0.9, width = 0.8) +
  # 在 total_rows 上方添加文本标签
  geom_text(aes(y = total_rows, label = paste(count_lt_adj, "/", total_rows, '\n', "(", sprintf("%.2f", proportion_adj * 100), "%)", sep = '')),
            vjust = -0.2, color = "black", size = 3) +
  # 设置主题和标签样式
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# 显示图形
p1


ggsave(filename = 'catalog_sig_proportion.pdf', plot = p1, width = 10, height = 6, dpi = 300)

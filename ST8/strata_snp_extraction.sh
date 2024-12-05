################################分层分析################################

#################################在两中心提取leadsnp#################################
R
rm(list=ls())
library(data.table)
library(dplyr)
#################################ZJ#################################
# 设置路径
input_dir <- "./Assoc_Lung_phenotype_low_fenceng/ZJ"
output_dir <- "./8.strata/1.data_extraction/ZJ"
snp_file <- "./8.strata/172_lead_snps_firstorder_shape.txt"

# 读取SNP位点信息
snps <- fread(snp_file)
# 创建组合变量
snps[, SNP_ID := paste(CHR, BP, sep = ":")]

# 获取所有以original开头且以.low.imputed.bolt.stats.gz为结尾的文件
file_list <- list.files(input_dir, pattern = "^original.*\\.low\\.imputed\\.bolt\\.stats\\.gz$", full.names = TRUE)

# 提取并保存数据
for (file in file_list) {
  data <- fread(file)
  data[, SNP_ID := paste(CHR, BP, sep = ":")]
  # 提取与SNP位点信息匹配的数据
  extracted_data <- data[SNP_ID %in% snps$SNP_ID]
  
  # 构建输出文件名
  base_name <- sub("\\.lung\\.low\\.imputed\\.bolt\\.stats\\.gz$", "", basename(file))
  output_file <- file.path(output_dir, paste0(base_name, "_172_lead_snps.txt"))

  fwrite(extracted_data, file = output_file, sep = "\t")
}


# 获取所有文件的列表
file_path <- "./8.strata/1.data_extraction/ZJ"
file_list <- list.files(path = file_path, pattern = "*.txt", full.names = TRUE)

# 函数：读取文件并添加phenotype列
read_and_add_phenotype <- function(file) {
  df <- fread(file)  # 使用fread读取文件
  phenotype <- gsub("_assoc.ZJ_172_lead_snps.txt$", "", basename(file))
  df[, phenotype := phenotype]
  return(df)
}

# 读取所有文件并合并到一个大数据框
data_list <- lapply(file_list, read_and_add_phenotype)
merged_df <- rbindlist(data_list)
fwrite(merged_df, file = "./8.strata/1.data_extraction/ZJ/merged_result_ZJ_all.csv")


# 按照SNP_ID去重，只保留P_BOLT_LMM最小的那一行，对于P值相同的重复行随机保留一行
result_df <- merged_df %>%
  group_by(SNP_ID) %>%
  filter(P_BOLT_LMM == min(P_BOLT_LMM)) %>%
  sample_n(1) %>%
  ungroup()

# 保存结果到新文件（可选）
fwrite(result_df, file = "./8.strata/1.data_extraction/ZJ/merged_result_ZJ.csv")
cat("合并和去重操作已完成。\n") #所有点在浙江的结果，性染色体的点较少直接grep看结果。




#################################JS#################################
# 设置路径
input_dir <- "./Assoc_Lung_phenotype_low_fenceng/JS"
output_dir <- "./8.strata/1.data_extraction/JS"
snp_file <- "./8.strata/172_lead_snps_firstorder_shape.txt"

# 读取SNP位点信息
snps <- fread(snp_file)
# 创建组合变量
snps[, SNP_ID := paste(CHR, BP, sep = ":")]

# 获取所有以original开头且以.low.imputed.bolt.stats.gz为结尾的文件
file_list <- list.files(input_dir, pattern = "^original.*\\.low\\.imputed\\.bolt\\.stats\\.gz$", full.names = TRUE)

# 提取并保存数据
for (file in file_list) {
  data <- fread(file)
  data[, SNP_ID := paste(CHR, BP, sep = ":")]
  extracted_data <- data[SNP_ID %in% snps$SNP_ID]
  base_name <- sub("\\.lung\\.low\\.imputed\\.bolt\\.stats\\.gz$", "", basename(file))
  output_file <- file.path(output_dir, paste0(base_name, "_172_lead_snps.txt"))
  fwrite(extracted_data, file = output_file, sep = "\t")
}



# 获取所有文件的列表
file_path <- "./8.strata/1.data_extraction/JS"
file_list <- list.files(path = file_path, pattern = "*.txt", full.names = TRUE)

# 函数：读取文件并添加phenotype列
read_and_add_phenotype <- function(file) {
  df <- fread(file)  # 使用fread读取文件
  phenotype <- gsub("_assoc.JS_172_lead_snps.txt$", "", basename(file))
  df[, phenotype := phenotype]
  return(df)
}

# 读取所有文件并合并到一个大数据框
data_list <- lapply(file_list, read_and_add_phenotype)
merged_df <- rbindlist(data_list)
fwrite(merged_df, file = "./8.strata/1.data_extraction/JS/merged_result_JS_all.csv") ###即为最终查找的数据框，根据位点和表型查找位点在ZJ人群的beta，se和p


# 按照SNP_ID去重，只保留P_BOLT_LMM最小的那一行，对于P值相同的重复行随机保留一行
result_df <- merged_df %>%
  group_by(SNP_ID) %>%
  filter(P_BOLT_LMM == min(P_BOLT_LMM)) %>%
  sample_n(1) %>%
  ungroup()

# 保存结果到新文件（可选）
fwrite(result_df, file = "./8.strata/1.data_extraction/JS/merged_result_JS.csv")

cat("合并和去重操作已完成。\n") #所有点在江苏的结果，性染色体的点较少直接grep看结果。





############################################吸烟非吸烟########################################
R
rm(list=ls())
library(data.table)
library(dplyr)
#################################非吸烟#################################
# 设置路径
input_dir <- "./Assoc_Lung_phenotype_low_fenceng/nonsmoke"
output_dir <- "./8.strata/1.data_extraction/nonsmoke"
snp_file <- "./8.strata/172_lead_snps_firstorder_shape.txt"

# 读取SNP位点信息
snps <- fread(snp_file)
# 创建组合变量
snps[, SNP_ID := paste(CHR, BP, sep = ":")]

# 获取所有以original开头且以.low.imputed.bolt.stats.gz为结尾的文件
file_list <- list.files(input_dir, pattern = "^original.*\\.low\\.imputed\\.bolt\\.stats\\.gz$", full.names = TRUE)

# 提取并保存数据
for (file in file_list) {
  # 读取文件
  data <- fread(file)
  # 创建组合变量
  data[, SNP_ID := paste(CHR, BP, sep = ":")]
  
  # 提取与SNP位点信息匹配的数据
  extracted_data <- data[SNP_ID %in% snps$SNP_ID]
  
  # 构建输出文件名
  base_name <- sub("\\.lung\\.low\\.imputed\\.bolt\\.stats\\.gz$", "", basename(file))
  output_file <- file.path(output_dir, paste0(base_name, "_172_lead_snps.txt"))
  
  # 保存提取的数据
  fwrite(extracted_data, file = output_file, sep = "\t")
}


# 文件路径
file_path <- "./8.strata/1.data_extraction/nonsmoke"

# 获取所有文件的列表
file_list <- list.files(path = file_path, pattern = "*.txt", full.names = TRUE)

# 函数：读取文件并添加phenotype列
read_and_add_phenotype <- function(file) {
  df <- fread(file)  # 使用fread读取文件
  phenotype <- gsub("_assoc.nonsmoke_172_lead_snps.txt$", "", basename(file))
  df[, phenotype := phenotype]
  return(df)
}

# 读取所有文件并合并到一个大数据框
data_list <- lapply(file_list, read_and_add_phenotype)
merged_df <- rbindlist(data_list)
fwrite(merged_df, file = "./8.strata/1.data_extraction/nonsmoke/merged_result_nonsmoke_all.csv")


# 按照SNP_ID去重，只保留P_BOLT_LMM最小的那一行，对于P值相同的重复行随机保留一行
result_df <- merged_df %>%
  group_by(SNP_ID) %>%
  filter(P_BOLT_LMM == min(P_BOLT_LMM)) %>%
  sample_n(1) %>%
  ungroup()

# 保存结果到新文件（可选）
fwrite(result_df, file = "./8.strata/1.data_extraction/nonsmoke/merged_result_nonsmoke.csv")

cat("合并和去重操作已完成。\n") #所有点在非吸烟的结果，性染色体的点较少直接grep看结果。


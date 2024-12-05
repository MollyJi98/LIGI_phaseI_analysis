#######################################计算膨胀系数(lambda1000)#######################################


R
rm(list = ls())

install.packages("writexl")
library(data.table)
library(writexl)

#p_value=phenotype_gwas$P_BOLT_LMM
#z = qnorm(p_value/ 2)
#lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)

file_list <- list.files("./Assoc_Lung_phenotype_low_3.8w_clean", full.names = TRUE)
# 创建一个空的数据框
result_df <- data.frame(file_name = character(0), lambda = numeric(0))
# 逐个读入文件并进行处理
for (file in file_list) {
  # 读入数据框，以文件名命名
  df <- fread(file)
  # 计算 z 值和 lambda
  p_value <- df$P_BOLT_LMM
  z <- qnorm(p_value / 2)
  lambda <- round(median(z^2, na.rm = TRUE) / 0.454, 3)
  # 将文件名和计算得到的lambda值添加到数据框中
  result_df <- rbind(result_df, data.frame(file_name = basename(file), lambda = lambda))
}

write_xlsx(result_df,"inflation.xlsx")






################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
################################对160个表型批量处理画QQ图################################
R
rm(list = ls())
library(CMplot)
library(data.table)




###############################################160个文件###############################################
### 获取文件列表
file_list <- list.files("./Assoc_Lung_phenotype_low_3.8w_clean", full.names = TRUE)
# 逐组读入文件并进行处理，5个一组用cmplot画qq图
for (i in seq(1, length(file_list), by=5)) {
  group_files <- file_list[i:(i+4)]
  # 创建一个空的数据框
  result_df <- data.frame(SNP = character(0) , CHR = character(0) , BP = character(0) )
  # 读入数据框并进行处理
  for (file in group_files) {
    # 读入数据框，以文件名命名
    df1 <- fread(file)
    # 保留第1, 2, 3, 16列，重命名第4列
    df2 <- df1[, c(1, 2, 3, 16)]
    colnames(df2)[4] <- gsub("original_|_assoc.result", "", basename(file))
    # 存储处理后的数据框
    result_df <- merge(result_df,df2,all=T)
  }
  # 生成文件名
  file_name <- gsub("original_shape_|original_firstorder_|_assoc.result", "", basename(group_files[1]))
  # 绘制QQ图
  CMplot(result_df, plot.type = "q", col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1","purple","grey"),
         multraits = TRUE, threshold = 1e-6, ylab.pos = 2, signal.pch = c(1,2,3,4,5), signal.cex = 1.2,
         signal.col = "red", conf.int = TRUE, box = FALSE, axis.cex = 1, file = "pdf", file.name = file_name,
         dpi = 300, file.output = TRUE, verbose = FALSE, width = 5, height = 5)
}


#####################后面在ai里拼图#####################
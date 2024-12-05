# IVW的结果整理
# 结果矩阵
rm(list = ls())

library(data.table)
library(tidyverse)

# 读入结果索引
files <- data.frame(
  names = system("ls /output_dir | grep mr.res.csv", intern = T)
) %>%
  mutate(
    res_path = paste0("/output_dir/", names)
  )

# 读入暴露索引
exp_list <- data.frame(
  exp_name = system("ls ./all_pheno_clumped", intern = T)
) %>%
  mutate(
    exp_name = str_extract(exp_name, ".*(?=\\.clumped\\.txt)")
  )

# 声明
res_com <- data.frame(
  exp_name = exp_list$exp_name
)
# loop
for (i in 1:nrow(files)) {
  res <- fread(files$res_path[i])
  out_col_name <- unique(res$outcome)
  snp_col_name <- paste0(out_col_name, "_nsnp")
  beta_col_name <- paste0(out_col_name, "_beta")
  p_col_name <- paste0(out_col_name, "_pval")
  se_col_name <- paste0(out_col_name, "_se")
  res <- res %>%
    filter(
      method == "Inverse variance weighted"
    ) %>%
    mutate(
      exp_name = str_extract(exposure, ".*(?=\\-)")
    ) %>%
    select(
      exp_name, nsnp, b, pval, se
    ) %>%
    rename(
      !!sym(snp_col_name) := nsnp,
      !!sym(beta_col_name) := b,
      !!sym(se_col_name) := se,
      !!sym(p_col_name) := pval
    )
  res_com <- merge(res_com, res, by = "exp_name", all = T)
  cat("-----", exp_list$exp_name[i], " has been processed ! ----- \n")
}
write.table(res_com, file = "res_matrix.txt", quote = F, row.names = F)

# 整理结果
rm(list = ls())

library(data.table)
library(tidyverse)

res <- fread("res_matrix.txt")
res <- as.data.frame(res)
df1 <- res %>%
  select(
    1, which(grepl("pval", colnames(res)))
  ) %>%
  pivot_longer(
    cols = ends_with("_pval"), # 选择所有表示 p 值的列
    names_to = "variable", # 新列名字
    values_to = "p_value" # 新列值
  ) %>%
  mutate(
    variable = str_replace(variable, "_pval", ""),
    merID = paste0(exp_name, ":", variable)
  )
df2 <- res %>%
  select(
    1, which(grepl("beta", colnames(res)))
  ) %>%
  pivot_longer(
    cols = ends_with("_beta"), # 选择所有表示 p 值的列
    names_to = "variable", # 新列名字
    values_to = "beta" # 新列值
  ) %>%
  mutate(
    variable = str_replace(variable, "_beta", ""),
    merID = paste0(exp_name, ":", variable)
  ) %>%
  select(
    merID, beta
  )
df3 <- res %>%
  select(
    1, which(grepl("_se", colnames(res)))
  ) %>%
  pivot_longer(
    cols = ends_with("_se"), # 选择所有表示 p 值的列
    names_to = "variable", # 新列名字
    values_to = "se" # 新列值
  ) %>%
  mutate(
    variable = str_replace(variable, "_se", ""),
    merID = paste0(exp_name, ":", variable)
  ) %>%
  select(
    merID, se
  )
df4 <- res %>%
  select(
    1, which(grepl("nsnp", colnames(res)))
  ) %>%
  pivot_longer(
    cols = ends_with("_nsnp"), # 选择所有表示 p 值的列
    names_to = "variable", # 新列名字
    values_to = "nsnp" # 新列值
  ) %>%
  mutate(
    variable = str_replace(variable, "_nsnp", ""),
    merID = paste0(exp_name, ":", variable)
  ) %>%
  select(
    merID, nsnp
  )

df <- merge(df1, df2, by = "merID")
df <- merge(df, df3, by = "merID")
df <- merge(df, df4, by = "merID")

df <- df %>%
  select(
    exp_name, variable, nsnp, beta, se, p_value
  )
library(readxl)
index <- read_xlsx("all_files_index.xlsx") 
index <- index[, c("name", "CN", "BBJ")]
index$BBJ <- str_extract(index$BBJ, "(?<=BBJ\\.).*(?=\\.v1)")
head(index)
df <- rename(df, BBJ = variable)
df <- merge(df, index, by = "BBJ", all.x = T)
df$p_fdr_adj <- p.adjust(df$p_value, method = "BH")
write.xlsx(df, file = "MR_Results_all.xlsx")
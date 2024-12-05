####################################同一表型两两分析遗传相关性####################################
# 设置路径
sumstats_dir="./step1_sumstats_result"
output_dir="./step5_phe_different_lobe_ldsc"
ldscores_dir="./ldsc/ldsc-master/eas_ldscores"

# 切换到目标文件夹
cd "$sumstats_dir"

# 遍历所有_munge.sumstats.gz文件，先排序再取5个一组
find . -name '*_munge.sumstats.gz' | sort | split -l 5 -d -a 3 - "$output_dir/input_files"

# 对每组文件夹进行遍历并创建文件夹
for folder in "$output_dir"/input_files*; do
    # 获取组内每个文件的完整路径，并将它们存储到数组中
    files=()
    while IFS= read -r file; do
        files+=("$file")
    done < "$folder"

    # 创建文件夹，并以第一个文件的 basename（去掉 "_CR_munge.sumstats.gz"）作为文件夹名
    first_file="${files[0]}"
    folder_name=$(basename "$first_file" | sed 's/_CR_munge.sumstats.gz//')
    mkdir -p "$output_dir/$folder_name"

    # 对于每组文件，计算两两之间的遗传相关性
    for ((i = 0; i < ${#files[@]}; i++)); do
        for ((j = i + 1; j < ${#files[@]}; j++)); do
            filename1=$(basename "${files[i]}" | sed 's/_munge.sumstats.gz//')
            filename2=$(basename "${files[j]}" | sed 's/_munge.sumstats.gz//')

            # 执行 ldsc.py 计算遗传相关性
            /home/jichen/ldsc-master/ldsc.py \
                --rg "${files[i]}","${files[j]}" \
                --ref-ld-chr "${ldscores_dir}/" \
                --w-ld-chr "${ldscores_dir}/" \
                --out "${output_dir}/${folder_name}/${filename1}.${filename2}_rg"
        done
    done
done



####################################整理结果####################################
# 加载必要的包
library(dplyr)

# Define the function to read the 62nd line of a file
read_62nd_line <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) >= 62) {
    return(lines[62])
  } else {
    return(NA)
  }
}

# Define the base directory
base_dir <- "./step5_phe_different_lobe_ldsc/"
sub_dirs <- list.dirs(base_dir, recursive = FALSE)
result_list <- list()

# Loop through each subdirectory
for (sub_dir in sub_dirs) {
  # Get list of all .log files in the current subdirectory
  log_files <- list.files(sub_dir, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)
  
  # Initialize an empty data frame to store the results for the current subdirectory
  sub_dir_data <- data.frame(line_62 = character(), filename = character(), stringsAsFactors = FALSE)
  
  # Loop through each log file and extract the 62nd line
  for (log_file in log_files) {
    line_62 <- read_62nd_line(log_file)
    filename <- sub("_rg\\.log$", "", basename(log_file))
    sub_dir_data <- rbind(sub_dir_data, data.frame(line_62 = line_62, filename = filename, stringsAsFactors = FALSE))
  }
  
  # Add the data frame for the current subdirectory to the result list
  result_list[[basename(sub_dir)]] <- sub_dir_data
}



library(dplyr)

# 定义读取文件第62行的函数（相关性结果在第62行）
read_62nd_line <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) >= 62) {
    return(lines[62])
  } else {
    return(NA)
  }
}

# 定义基本目录
base_dir <- "./step5_phe_different_lobe_ldsc/"

# 获取基本目录下所有子目录的列表
sub_dirs <- list.dirs(base_dir, recursive = FALSE)
sub_dirs <- sub_dirs[grepl("^original_", basename(sub_dirs))]

# 初始化一个空的列表来存储每个子目录的结果
result_list <- list()

# 遍历每个子目录
for (sub_dir in sub_dirs) {
  # 获取当前子目录下所有的.log文件
  log_files <- list.files(sub_dir, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)
  
  # 初始化一个空的数据框来存储当前子目录的结果
  sub_dir_data <- data.frame(line_62 = character(), filename = character(), stringsAsFactors = FALSE)
  
  # 遍历每个.log文件并提取第62行
  for (log_file in log_files) {
    line_62 <- read_62nd_line(log_file)
    filename <- sub("_rg\\.log$", "", basename(log_file))
    sub_dir_data <- rbind(sub_dir_data, data.frame(line_62 = line_62, filename = filename, stringsAsFactors = FALSE))
  }
  
  # 将当前子目录的数据框添加到结果列表中
  result_list[[basename(sub_dir)]] <- sub_dir_data
}

# 处理结果列表中的每个数据框
for (name in names(result_list)) {
  df <- result_list[[name]]
  
  # 使用基础R的strsplit函数按空格分割每一行并创建新的列
  split_lines <- strsplit(df$line_62, "\\s+")
  split_df <- do.call(rbind, lapply(split_lines, function(x) {
    if (length(x) == 12) {
      return(x)
    } else {
      return(rep(NA, 12))  # 如果元素个数不足12个，则用NA填充
    }
  }))
  
  colnames(split_df) <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")
  df <- cbind(df, split_df)
  df <- df %>% select(-line_62)  # 删除原始的line_62列
  
  # 更新结果列表中的数据框
  result_list[[name]] <- df
}


#########################################画热图#########################################å
# 加载必要的包
library(ggplot2)
library(gridExtra)

# 定义绘制热图的函数
plot_heatmap <- function(corr_matrix, pval_matrix, title) {
  # 去掉行名和列名中的 `_` 及其前面的内容
  rownames(corr_matrix) <- gsub(".*_", "", rownames(corr_matrix))
  colnames(corr_matrix) <- gsub(".*_", "", colnames(corr_matrix))
  rownames(pval_matrix) <- gsub(".*_", "", rownames(pval_matrix))
  colnames(pval_matrix) <- gsub(".*_", "", colnames(pval_matrix))
  
  # 将矩阵转换为长格式
  melted_corr_matrix <- data.frame(
    Phenotype1 = rep(rownames(corr_matrix), each = ncol(corr_matrix)),
    Phenotype2 = rep(colnames(corr_matrix), times = nrow(corr_matrix)),
    Correlation = as.vector(corr_matrix),
    PValue = as.vector(pval_matrix)
  )

  # 创建标注信息
  melted_corr_matrix$label <- ifelse(is.na(melted_corr_matrix$Correlation), "NA",
                                     ifelse(melted_corr_matrix$PValue < 0.05, "*", ""))

  melted_corr_matrix$Correlation[is.na(melted_corr_matrix$Correlation)] <- NA

  ggplot(melted_corr_matrix, aes(x = Phenotype1, y = Phenotype2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(round(Correlation, 2), ifelse(label == "NA", "", label))),
              color = "black", size = 3) +
    scale_fill_gradient2(low = "#5b98e0", high = "firebrick3", mid = "white", midpoint = 0,
                         limit = c(-1, 1), space = "Lab", name = "rg") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),  # 隐藏 x 轴标题
      axis.title.y = element_blank(),  # 隐藏 y 轴标题
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 15, face = "bold"),  # 设置标题字体大小
      axis.text = element_text(size = 15),   # 设置轴文本字体大小
      legend.text = element_text(size = 10),  # 设置图例文本字体大小
      legend.title = element_text(size = 10)  # 设置图例标题字体大小
    ) +
    ggtitle(title)
}

# 初始化一个空列表来存储热图对象
heatmaps <- list()

# 遍历 result_list 中的每个数据框
for (name in names(result_list)) {
  df <- result_list[[name]]

  # 去掉 p1 和 p2 列中的 ./original_、_munge.sumstats.gz、firstorder_ 以及 shape_
  df$p1 <- gsub("./original_|_munge.sumstats.gz|firstorder_|shape_", "", df$p1)
  df$p2 <- gsub("./original_|_munge.sumstats.gz|firstorder_|shape_", "", df$p2)

  # 替换 RobustMeanAbsoluteDeviation 为 RMAD，MeanAbsoluteDeviation 为 MAD
  df$p1 <- gsub("RobustMeanAbsoluteDeviation", "RMAD", df$p1)
  df$p2 <- gsub("RobustMeanAbsoluteDeviation", "RMAD", df$p2)
  df$p1 <- gsub("MeanAbsoluteDeviation", "MAD", df$p1)
  df$p2 <- gsub("MeanAbsoluteDeviation", "MAD", df$p2)

  # 过滤掉 p1 和 p2 为 NA 的行
  df <- df[!is.na(df$p1) & !is.na(df$p2),]

  # 创建一个空的相关性矩阵和 p 值矩阵
  phenotypes <- unique(c(df$p1, df$p2))
  corr_matrix <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
  pval_matrix <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
  colnames(corr_matrix) <- phenotypes
  rownames(corr_matrix) <- phenotypes
  colnames(pval_matrix) <- phenotypes
  rownames(pval_matrix) <- phenotypes

  # 填充相关性矩阵和 p 值矩阵
  for (i in 1:nrow(df)) {
    p1 <- df$p1[i]
    p2 <- df$p2[i]
    rg <- as.numeric(df$rg[i])
    p <- as.numeric(df$p[i])
    corr_matrix[p1, p2] <- ifelse(is.na(rg), 0, ifelse(rg < -1, -1, ifelse(rg > 1, 1, rg)))
    corr_matrix[p2, p1] <- ifelse(is.na(rg), 0, ifelse(rg < -1, -1, ifelse(rg > 1, 1, rg)))  # 矩阵对称填充
    pval_matrix[p1, p2] <- ifelse(is.na(p), NA, p)
    pval_matrix[p2, p1] <- ifelse(is.na(p), NA, p)  # 矩阵对称填充
  }

  # 将对角线上的值设为1（自身相关性）
  diag(corr_matrix) <- 1
  diag(pval_matrix) <- 0  # 对角线上的 p 值设为 0

  # 绘制相关性矩阵热图并添加到热图列表中
  plot_title <- gsub("Correlation Heatmap for ", "", gsub("original_", "", name))
  heatmap_plot <- plot_heatmap(corr_matrix, pval_matrix, plot_title)

  # 保存单个热图
  ggsave(paste0("./step5_phe_different_lobe_ldsc/plot/", plot_title, "_correlation_heatmap.pdf"),
         plot = heatmap_plot, width = 5, height = 5, device = cairo_pdf)

  # 添加热图到列表
  heatmaps[[name]] <- heatmap_plot
}

# 定义每页包含的行列数
rows_per_page <- 4
cols_per_page <- 4

# 计算页数
num_plots <- length(heatmaps)
plots_per_page <- rows_per_page * cols_per_page
num_pages <- ceiling(num_plots / plots_per_page)

# 按页保存热图
for (i in 1:num_pages) {
  start_index <- (i - 1) * plots_per_page + 1
  end_index <- min(i * plots_per_page, num_plots)
  current_page_plots <- heatmaps[start_index:end_index]

  # 将当前页的热图拼接成一个 4x4 的图网格
  current_combined_plot <- marrangeGrob(grobs = current_page_plots, ncol = cols_per_page, nrow = rows_per_page, top = NULL)

  # 保存当前页的拼接热图
  ggsave(paste0("./step5_phe_different_lobe_ldsc/plot/combined_correlation_heatmap_page_", i, ".pdf"),
         plot = current_combined_plot, width = 16, height = 16, device = cairo_pdf, limitsize = FALSE)
}

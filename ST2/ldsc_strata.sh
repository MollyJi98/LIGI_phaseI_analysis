########################################20240513 ldsc for strata analysis by JC########################################

########################################step 1. 针对所有表型保留质控后位点，配rs号########################################
cd ./
mkdir ./ldsc #创建ldsc路径
mkdir ./ldsc/ldsc_step1 #创建ldsc step1路径
mkdir ./ldsc/ldsc_step1/ZJ
mkdir ./ldsc/ldsc_step1/JS
mkdir ./ldsc/ldsc_step1/nonsmoke
mkdir ./ldsc/ldsc_step1/smoke

##注意，以下所有r代码写到r的脚本里，用.r命名后，用dos2unix转换成unix脚本后，运行再用sh脚本运行##


###################ZJ 数据质控及配rs号（step1_ZJ.r）###################
R
library(data.table)
setwd('./ZJ')

# 设置文件路径和文件名模式
file_path <- "./ZJ" #输入文件路径
ldsc_step1_path <- "./ldsc/ldsc_step1/ZJ/" #输出文件路径
file_pattern <- "ZJ.lung.low.imputed.bolt.stats.gz" #输入文件指定后缀

# 读取 SNP_chrpos_matching.txt 文件，用于匹配
link <- fread('./ldsc/SNP_chrpos_matching.txt')

# 获取所有以 ZJ.lung.low.imputed.bolt.stats.gz 为结尾的文件
files <- list.files(file_path, pattern = file_pattern, full.names = TRUE)

# 循环处理每个文件
for (file in files) {
  # 读取文件
  a <- fread(file)
  # 生成 chrpos 列
  a$chrpos <- paste(a$CHR, a$BP, sep = ":")
  # 保留与 link 中的 chrpos 匹配的行
  a1 <- subset(a, a$chrpos %in% link$chrpos)
  # 合并 rsid
  a2 <- merge(a1, link, by = "chrpos")
  
  # 调整列的顺序
  a3 <- a2[, c("SNP", "CHR", "BP", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "CHISQ_LINREG", "P_LINREG", 
               "BETA", "SE", "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "P_BOLT_LMM", "chrpos", "RSID")]
  
  # 生成输出文件名
  output_file_name <- paste0(ldsc_step1_path, gsub(".ZJ.lung.low.imputed.bolt.stats.gz", "_for_LDSC.txt", basename(file)))
  
  # 写入文件
  write.table(a3, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
}

##!/bin/sh
##cd ./ldsc/ldsc_step1/JS
#module load bioinfo
#module load R
#Rscript step1_JS.r
###################JS 数据质控及配rs号（step1_JS.r）###################
R
library(data.table)
setwd('./JS')

# 设置文件路径和文件名模式
file_path <- "./JS" #输入文件路径
ldsc_step1_path <- "./ldsc/ldsc_step1/JS/" #输出文件路径
file_pattern <- "JS.lung.low.imputed.bolt.stats.gz" #输入文件指定后缀

# 读取 SNP_chrpos_matching.txt 文件，用于匹配
link <- fread('./ldsc/SNP_chrpos_matching.txt')

# 获取所有以 JS.lung.low.imputed.bolt.stats.gz 为结尾的文件
files <- list.files(file_path, pattern = file_pattern, full.names = TRUE)

# 循环处理每个文件
for (file in files) {
  # 读取文件
  a <- fread(file)
  
  # 生成 chrpos 列
  a$chrpos <- paste(a$CHR, a$BP, sep = ":")
  
  # 保留与 link 中的 chrpos 匹配的行
  a1 <- subset(a, a$chrpos %in% link$chrpos)
  
  # 合并 rsid
  a2 <- merge(a1, link, by = "chrpos")
  
  # 调整列的顺序
  a3 <- a2[, c("SNP", "CHR", "BP", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "CHISQ_LINREG", "P_LINREG", 
               "BETA", "SE", "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "P_BOLT_LMM", "chrpos", "RSID")]
  
  # 生成输出文件名
  output_file_name <- paste0(ldsc_step1_path, gsub(".JS.lung.low.imputed.bolt.stats.gz", "_for_LDSC.txt", basename(file)))
  
  # 写入文件
  write.table(a3, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
}



###################smoke 数据质控及配rs号(step1_smoke.r)###################
R
library(data.table)
setwd('./smoke')

# 设置文件路径和文件名模式
file_path <- "./smoke" #输入文件路径
ldsc_step1_path <- "./ldsc/ldsc_step1/smoke/" #输出文件路径
file_pattern <- "smoke.lung.low.imputed.bolt.stats.gz" #输入文件指定后缀

# 读取 SNP_chrpos_matching.txt 文件，用于匹配
link <- fread('./ldsc/SNP_chrpos_matching.txt')

# 获取所有以 smoke.lung.low.imputed.bolt.stats.gz 为结尾的文件
files <- list.files(file_path, pattern = file_pattern, full.names = TRUE)

# 循环处理每个文件
for (file in files) {
  # 读取文件
  a <- fread(file)
  
  # 生成 chrpos 列
  a$chrpos <- paste(a$CHR, a$BP, sep = ":")
  
  # 保留与 link 中的 chrpos 匹配的行
  a1 <- subset(a, a$chrpos %in% link$chrpos)
  
  # 合并 rsid
  a2 <- merge(a1, link, by = "chrpos")
  
  # 调整列的顺序
  a3 <- a2[, c("SNP", "CHR", "BP", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "CHISQ_LINREG", "P_LINREG", 
               "BETA", "SE", "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "P_BOLT_LMM", "chrpos", "RSID")]
  
  # 生成输出文件名
  output_file_name <- paste0(ldsc_step1_path, gsub(".smoke.lung.low.imputed.bolt.stats.gz", "_for_LDSC.txt", basename(file)))
  
  # 写入文件
  write.table(a3, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
}


###################nonsmoke 数据质控及配rs号(step1_nonsmoke.r)###################
R
library(data.table)
setwd('./nonsmoke')

# 设置文件路径和文件名模式
file_path <- "./nonsmoke" #输入文件路径
ldsc_step1_path <- "./ldsc/ldsc_step1/nonsmoke/" #输出文件路径
file_pattern <- "nonsmoke.lung.low.imputed.bolt.stats.gz" #输入文件指定后缀

# 读取 SNP_chrpos_matching.txt 文件，用于匹配
link <- fread('./ldsc/SNP_chrpos_matching.txt')

# 获取所有以 nonsmoke.lung.low.imputed.bolt.stats.gz 为结尾的文件
files <- list.files(file_path, pattern = file_pattern, full.names = TRUE)

# 循环处理每个文件
for (file in files) {
  # 读取文件
  a <- fread(file)
  
  # 生成 chrpos 列
  a$chrpos <- paste(a$CHR, a$BP, sep = ":")
  
  # 保留与 link 中的 chrpos 匹配的行
  a1 <- subset(a, a$chrpos %in% link$chrpos)
  
  # 合并 rsid
  a2 <- merge(a1, link, by = "chrpos")
  
  # 调整列的顺序
  a3 <- a2[, c("SNP", "CHR", "BP", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "CHISQ_LINREG", "P_LINREG", 
               "BETA", "SE", "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "P_BOLT_LMM", "chrpos", "RSID")]
  
  # 生成输出文件名
  output_file_name <- paste0(ldsc_step1_path, gsub(".nonsmoke.lung.low.imputed.bolt.stats.gz", "_for_LDSC.txt", basename(file)))
  
  # 写入文件
  write.table(a3, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
}







########################################step 2. 运行munge脚本进行数据整理########################################
mkdir ./ldsc #创建ldsc路径
mkdir ./ldsc/ldsc_step2 #创建ldsc step2路径
mkdir ./ldsc/ldsc_step2/ZJ
mkdir ./ldsc/ldsc_step2/JS
mkdir ./ldsc/ldsc_step2/nonsmoke
mkdir ./ldsc/ldsc_step2/smoke


#用munge_sumstats进行处理#
source activate ldsc

##############################批量写循环（ZJ）##############################
#cd 到输出文件的路径
cd ./ldsc/ldsc_step2/ZJ
# 设置输入文件路径
input_dir="./ldsc/ldsc_step1/ZJ"

# 获取文件列表，以original_开头，_for_LDSC.txt为结尾的文件
file_list=$(ls "$input_dir"/original_*_for_LDSC.txt)

# 循环处理每个文件
for file in $file_list; do
    # 提取文件名（去除路径和文件扩展名）
    filename=$(basename "$file" _for_LDSC.txt)

    # 提取指定列并重命名列名，输出到临时文件
    awk -F'\t' 'BEGIN{OFS="\t"; print "CHR", "BP", "A1", "A2", "BETA", "P", "FRQ", "INFO", "SNP"} NR>1 {printf "%d\t%d\t%s\t%s\t%.6f\t%.200f\t%.3f\t%.3f\t%s\n", $2, $3, $5, $6, $11+0, $16+0, $7+0, $8+0, $18}' "$file" \
    > output.txt

    # 设置输出文件名
    output_name="${filename}_munge"

    # 运行 munge_sumstats.py 命令(！！！！！注意这里要指定为实际的ldsc的py脚本路径！！！！！)
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/munge_sumstats.py --N 25100 --info-min 0.3 --sumstats output.txt --out "${output_name}"

    # 删除临时文件
    rm output.txt
done


##############################批量写循环（JS）##############################
#cd 到输出文件的路径
cd ./ldsc/ldsc_step2/JS
# 设置输入文件路径
input_dir="./ldsc/ldsc_step1/JS"

# 获取文件列表，以original_开头，_for_LDSC.txt为结尾的文件
file_list=$(ls "$input_dir"/original_*_for_LDSC.txt)

# 循环处理每个文件
for file in $file_list; do
    # 提取文件名（去除路径和文件扩展名）
    filename=$(basename "$file" _for_LDSC.txt)

    # 提取指定列并重命名列名，输出到临时文件
    awk -F'\t' 'BEGIN{OFS="\t"; print "CHR", "BP", "A1", "A2", "BETA", "P", "FRQ", "INFO", "SNP"} NR>1 {printf "%d\t%d\t%s\t%s\t%.6f\t%.200f\t%.3f\t%.3f\t%s\n", $2, $3, $5, $6, $11+0, $16+0, $7+0, $8+0, $18}' "$file" \
    > output.txt

    # 设置输出文件名
    output_name="${filename}_munge"

    # 运行 munge_sumstats.py 命令(！！！！！注意这里要指定为实际的ldsc的py脚本路径！！！！！)
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/munge_sumstats.py --N 10369 --info-min 0.3 --sumstats output.txt --out "${output_name}"

    # 删除临时文件
    rm output.txt
done


##############################批量写循环（smoke）##############################
#cd 到输出文件的路径
cd ./ldsc/ldsc_step2/smoke
# 设置输入文件路径
input_dir="./ldsc/ldsc_step1/smoke"

# 获取文件列表，以original_开头，_for_LDSC.txt为结尾的文件
file_list=$(ls "$input_dir"/original_*_for_LDSC.txt)

# 循环处理每个文件
for file in $file_list; do
    # 提取文件名（去除路径和文件扩展名）
    filename=$(basename "$file" _for_LDSC.txt)

    # 提取指定列并重命名列名，输出到临时文件
    awk -F'\t' 'BEGIN{OFS="\t"; print "CHR", "BP", "A1", "A2", "BETA", "P", "FRQ", "INFO", "SNP"} NR>1 {printf "%d\t%d\t%s\t%s\t%.6f\t%.200f\t%.3f\t%.3f\t%s\n", $2, $3, $5, $6, $11+0, $16+0, $7+0, $8+0, $18}' "$file" \
    > output.txt

    # 设置输出文件名
    output_name="${filename}_munge"

    # 运行 munge_sumstats.py 命令(！！！！！注意这里要指定为实际的ldsc的py脚本路径！！！！！)
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/munge_sumstats.py --N 15313 --info-min 0.3 --sumstats output.txt --out "${output_name}"

    # 删除临时文件
    rm output.txt
done


##############################批量写循环（nonsmoke）##############################
#cd 到输出文件的路径
cd ./ldsc/ldsc_step2/nonsmoke
# 设置输入文件路径
input_dir="./ldsc/ldsc_step1/nonsmoke"

# 获取文件列表，以original_开头，_for_LDSC.txt为结尾的文件
file_list=$(ls "$input_dir"/original_*_for_LDSC.txt)

# 循环处理每个文件
for file in $file_list; do
    # 提取文件名（去除路径和文件扩展名）
    filename=$(basename "$file" _for_LDSC.txt)

    # 提取指定列并重命名列名，输出到临时文件
    awk -F'\t' 'BEGIN{OFS="\t"; print "CHR", "BP", "A1", "A2", "BETA", "P", "FRQ", "INFO", "SNP"} NR>1 {printf "%d\t%d\t%s\t%s\t%.6f\t%.200f\t%.3f\t%.3f\t%s\n", $2, $3, $5, $6, $11+0, $16+0, $7+0, $8+0, $18}' "$file" \
    > output.txt

    # 设置输出文件名
    output_name="${filename}_munge"

    # 运行 munge_sumstats.py 命令(！！！！！注意这里要指定为实际的ldsc的py脚本路径！！！！！)
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/munge_sumstats.py --N 20148 --info-min 0.3 --sumstats output.txt --out "${output_name}"

    # 删除临时文件
    rm output.txt
done



########################################step 3. 计算遗传相关性########################################
conda activate ldsc
mkdir ./ldsc/ldsc_step3/ZJ_JS
mkdir ./ldsc/ldsc_step3/smoke_nonsmoke


###########ZJ&JS#################
#!/bin/bash
source activate ldsc

cd ./ldsc/ldsc_step3/ZJ_JS

ldsc_step2_path="./ldsc/ldsc_step2"
file_a_path="$ldsc_step2_path/ZJ"
file_b_path="$ldsc_step2_path/JS"
output_dir="$(pwd)"
ldscores_dir="./Software/ldsc-master/eas_ldscores"

file_a=$(mktemp)
file_b=$(mktemp)

# 将文件列表存储到临时文件中
find "$file_a_path" -maxdepth 1 -type f -name '*.gz' | sort > "$file_a"
find "$file_b_path" -maxdepth 1 -type f -name '*.gz' | sort > "$file_b"

# 合并两个文件的内容并输出到新文件
paste "$file_a" "$file_b" > merged_files.txt

# 逐行遍历合并后的文件，赋值给 filename1 和 filename2
while IFS= read -r line; do
    filename1=$(echo "$line" | awk '{print $1}')
    filename2=$(echo "$line" | awk '{print $2}')
    a="${filename1##*/}"
    a="${a%%_munge.sumstats.gz}"
    b="${filename2##*/}"
    b="${b%%_munge.sumstats.gz}"
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/ldsc.py \
    --rg "$filename1","$filename2" \
    --ref-ld-chr "${ldscores_dir}/" \
    --w-ld-chr "${ldscores_dir}/" \
    --out "${output_dir}/${a}.${b}_rg"
done < merged_files.txt

# 删除临时文件
rm "$file_a" "$file_b" merged_files.txt


###########smoke&nonsmoke#################
#!/bin/bash

source activate ldsc

cd ./ldsc/ldsc_step3/smoke_nonsmoke

ldsc_step2_path="./ldsc/ldsc_step2"
file_a_path="$ldsc_step2_path/smoke"
file_b_path="$ldsc_step2_path/nonsmoke"
output_dir="$(pwd)"
ldscores_dir="./Software/ldsc-master/eas_ldscores"

file_a=$(mktemp)
file_b=$(mktemp)

# 将文件列表存储到临时文件中
find "$file_a_path" -maxdepth 1 -type f -name '*.gz' | sort > "$file_a"
find "$file_b_path" -maxdepth 1 -type f -name '*.gz' | sort > "$file_b"

# 合并两个文件的内容并输出到新文件
paste "$file_a" "$file_b" > merged_files.txt

# 逐行遍历合并后的文件，赋值给 filename1 和 filename2
while IFS= read -r line; do
    filename1=$(echo "$line" | awk '{print $1}')
    filename2=$(echo "$line" | awk '{print $2}')
    a="${filename1##*/}"
    a="${a%%_munge.sumstats.gz}"
    b="${filename2##*/}"
    b="${b%%_munge.sumstats.gz}"
    ./software/anaconda3/envs/ldsc/bin/python ./Software/ldsc-master/ldsc.py \
    --rg "$filename1","$filename2" \
    --ref-ld-chr "${ldscores_dir}/" \
    --w-ld-chr "${ldscores_dir}/" \
    --out "${output_dir}/${a}.${b}_rg"
done < merged_files.txt

# 删除临时文件
rm "$file_a" "$file_b" merged_files.txt




##################################跑完的所有数据储存在硬盘，整理结果##################################
# 加载必要的库
library(dplyr)
library(stringr)

# 定义数据目录
input_dir <- './ldsc_step3/smoke_nonsmoke'

# 获取所有log文件路径
file_paths <- list.files(input_dir, pattern = "\\.log$", full.names = TRUE)

# 初始化一个空的数据框用于存储结果
results_df <- data.frame(file_name = character(), genetic_correlation = character(), p_value = character(), stringsAsFactors = FALSE)

# 遍历每个文件进行处理
for (file_path in file_paths) {
  # 读取文件内容
  file_content <- readLines(file_path)
  
  # 查找包含"Genetic Correlation: "的行
  gc_line <- file_content[grep("Genetic Correlation: ", file_content)]
  
  # 查找包含"P: "的行
  p_line <- file_content[grep("P: ", file_content)]
  
  # 提取"Genetic Correlation: "到" ("之间的内容
  gc_value <- if (length(gc_line) > 0) {
    str_extract(gc_line, "(?<=Genetic Correlation: ).*?(?= \\()")
  } else {
    NA
  }
  
  # 提取"P: "之后的内容
  p_value <- if (length(p_line) > 0) {
    str_extract(p_line, "(?<=P: ).*")
  } else {
    NA
  }
  
  # 提取文件名
  file_name <- basename(file_path)
  
  # 将结果添加到数据框中
  results_df <- results_df %>% 
    add_row(file_name = file_name, genetic_correlation = gc_value, p_value = p_value)
}

# 将结果保存到CSV文件
output_file <- 'genetic_correlations_smoke_nonsmoke.csv'
write.csv(results_df, output_file, row.names = FALSE)



# 定义数据目录
input_dir <- './ldsc_step3/ZJ_JS'

# 获取所有log文件路径
file_paths <- list.files(input_dir, pattern = "\\.log$", full.names = TRUE)

# 初始化一个空的数据框用于存储结果
results_df <- data.frame(file_name = character(), genetic_correlation = character(), p_value = character(), stringsAsFactors = FALSE)

# 遍历每个文件进行处理
for (file_path in file_paths) {
  # 读取文件内容
  file_content <- readLines(file_path)
  
  # 查找包含"Genetic Correlation: "的行
  gc_line <- file_content[grep("Genetic Correlation: ", file_content)]
  
  # 查找包含"P: "的行
  p_line <- file_content[grep("P: ", file_content)]
  
  # 提取"Genetic Correlation: "到" ("之间的内容
  gc_value <- if (length(gc_line) > 0) {
    str_extract(gc_line, "(?<=Genetic Correlation: ).*?(?= \\()")
  } else {
    NA
  }
  
  # 提取"P: "之后的内容
  p_value <- if (length(p_line) > 0) {
    str_extract(p_line, "(?<=P: ).*")
  } else {
    NA
  }
  
  # 提取文件名
  file_name <- basename(file_path)
  
  # 将结果添加到数据框中
  results_df <- results_df %>% 
    add_row(file_name = file_name, genetic_correlation = gc_value, p_value = p_value)
}

# 将结果保存到CSV文件
output_file <- 'genetic_correlations_ZJ_JS.csv'
write.csv(results_df, output_file, row.names = FALSE)

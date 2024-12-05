########################clumping for lead snp鉴定每个表型的lead snp################

###################################################step1##########################################
# 设置输入和输出路径
input_dir="./0_sig_remove_outlier_snps/5e-8"  #gwas significcant data(<5e-8)
output_dir="./lead_snps_160phe" #desired output dir

# 遍历输入路径下所有的txt文件
for file in "$input_dir"/*.txt; do
    # 提取文件名（不含路径和扩展名）
    filename=$(basename "$file" .txt)
    # 提取第18和16列，写入临时文件并修改列名
    awk 'BEGIN {OFS="\t"} {if(NR==1) print "SNP", "P"; else print $18, $16}' "$file" > "${filename}_temp.txt"
    # 构造输出文件名
    output_file="${filename/_5e-8}_clumping_0.1"
    # 执行 clumping
    plink --bfile ./chr_all_1kgv3_2015_EAS \
          --clump "${filename}_temp.txt" \
          --clump-kb 1000 \
          --clump-p1 0.05 \
          --clump-p2 0.05 \
          --clump-r2 0.1 \
          --out "${output_dir}/${output_file}"
done




####################################step2(table4&6)#############################################
R
# 设置输入和输出路径
input_dir <- "./lead_snps_160phe"
output_file <- "lead_snps_160_phenotypes.txt"

# 获取目录中所有.clumped为后缀的文件列表
file_list <- list.files(input_dir, pattern = "\\.clumped$", full.names = TRUE)

# 初始化一个空数据框来存储合并后的结果
merged_data <- data.frame()

# 循环读取每个文件并合并到数据框中
for (file in file_list) {
  # 读取文件
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  
  # 提取文件名并去除路径和文件扩展名
  file_name <- basename(file)
  file_name <- gsub("_assoc_sig_clumping_0\\.1\\.clumped$", "", file_name)
  
  # 添加文件名列
  data$File_Name <- file_name
  
  # 合并到结果数据框中
  merged_data <- rbind(merged_data, data)
}

# 将结果保存到文件中
write.table(merged_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


#annotion
perl ./annovar/table_annovar.pl \
    lead_snps_160phe_for_annovar.txt \
    ./annovar/humandb/ \
    -buildver hg19 \
    -out lead_snps_160phe_for_annovar_res \
    -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -csvout -polish








####################################step3 reclumping(table5&7)#############################################
#first order
plink --bfile ./chr_all_1kgv3_2015_EAS \
      --clump lead_snp_firstorder_for_clumping2.txt \
      --clump-kb 1000 \
      --clump-p1 0.05 \
      --clump-r2 0.1 \
      --out lead_snp_firstorder_for_clumping2_res


#shape
plink --bfile ./chr_all_1kgv3_2015_EAS \
      --clump lead_snp_shape_for_clumping2.txt \
      --clump-kb 1000 \
      --clump-p1 0.05 \
      --clump-r2 0.1 \
      --out lead_snp_shape_for_clumping2_res
 

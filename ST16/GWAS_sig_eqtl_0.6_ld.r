########################900lead snp结果提取高ld(>0.6)位点########################
# 取所有lead snp以及其高ld(0.6)的位点
# lead snp共390个

#!/bin/bash

input_file="390_lead_snps.txt"
output_dir="./ld_0.6_extract"
bfile="./chr_all_1kgv3_2015_EAS"


while IFS= read -r rs_id; do
    # 去除行末空格
    rs_id=$(echo "$rs_id" | tr -d '[:space:]')
    
    # 检查是否是空行
    if [ -z "$rs_id" ]; then
        continue
    fi

    # 执行 PLINK 命令
    plink --bfile "$bfile" \
          --r2 \
          --ld-snp "$rs_id" \
          --ld-window-kb 1000 \
          --ld-window 99999 \
          --ld-window-r2 0.6 \
          --out "${output_dir}/${rs_id}_LD_snps"
done < "$input_file"

#清点数量
ls -l ./ld_0.6_extract/*_snps.ld | wc -l #388

#把所有文件合并到一个文件里
cd ./ld_0.6_extract
find ./ld_0.6_extract -name "*.ld" -exec cat {} + > 388_leadsnps_merged_data_ld_0.6.txt 


R
rm(list=ls())
setwd('./ld_0.6_extract')

library(data.table)
library(dplyr)
library(purrr)

#读取ld关系和eqtl数据
ld_data<-fread('388_leadsnps_merged_data_ld_0.6.txt')
eqtl<-fread('eqtl_RES_unique_snpNgene_withMHC.txt')

#生成chr:bp列方便匹配
ld_data$chr_bp<-paste(ld_data$"CHR_B",ld_data$"BP_B",sep=":")
eqtl$chr_bp<-paste(eqtl$"CHR.x",eqtl$BP,sep=":")

#挑选的ld关系数据里的位点
eqtl2<-subset(eqtl,eqtl$chr_bp %in% ld_data$chr_bp) #91645
#取出fdr矫正显著的eqtl
eqtl3<-subset(eqtl2,eqtl2$qvalue<0.05)
write.table(eqtl3,'388_ld_0.6_eqtl_RES_unique_snpNgene_sig_withMHC.txt', quote=F ,sep = "\t", row.names = FALSE)


#对于leadsnp的所有位点，按照每个phenotype进行处理，对lead snp合并高ld位点，最终得到每个phenotype应该包含的位点
#即最终包含的是每个表型lead snp及其高ld位点
lead_snp<-fread('908-3_lead_snps.txt')

#对于lead_snp数据框的每一种phenotype分别处理，对于每一种phenotype，对于依照SNP这一列与ld_data的SNP_A这一列进行多对多的合并，最终得到每一种phenotype的leadsnp以及与其有ld的关系的snp的数据框
result <- lead_snp %>%
  split(.$phenotype) %>%
  map_dfr(~ {
    phenotype_data <- .
    phenotype_ld <- ld_data %>%
      inner_join(phenotype_data, by = c("SNP_A" = "SNP"))
    phenotype_ld$phenotype <- unique(phenotype_data$phenotype)
    phenotype_ld
  })

write.table(result,'908-3_ld_0.6_lead_snps_for_merge.txt', quote=F ,sep = "\t", row.names = FALSE)
# !!以上得到的res就是每一种phenotype的leadsnp以及与其有ld的关系的snp的数据框，但不确保这些点都显著！！

#对于result的每一种phenotype分别处理，提取每种phenotype的行，在gwas_sig数据框中找到phenotype和File_Name一致的行，并对每种phenotype的行按照chrpos识别，只保留在gwas_sig对应表型中也出现的chrpos
gwas_sig<-fread('all_traits_sig_merge_new.txt') #读取gwas显著的数据框


result_filtered <- result %>%
  group_by(phenotype) %>%
  do({
    # 提取当前 phenotype 的数据
    current_result <- .
    
    # 在 gwas_sig 中找到对应的 File_Name
    file_name <- gwas_sig %>% 
      filter(File_Name == current_result$phenotype[1]) %>% 
      pull(File_Name) %>% 
      unique()
    
    # 如果找到了匹配的 File_Name
    if (length(file_name) > 0) {
      # 在 gwas_sig 中找到对应的 chrpos
      valid_chrpos <- gwas_sig %>%
        filter(File_Name == current_result$phenotype[1]) %>%
        select(chrpos,ALLELE1,ALLELE0, BETA, SE, P_BOLT_LMM)
      
      # 只保留在 gwas_sig 中也出现的 chrpos，并合并 BETA, SE, 和 P_BOLT_LMM 列
      filtered_result <- current_result %>%
        filter(chr_bp %in% valid_chrpos$chrpos) %>%
        left_join(valid_chrpos, by = c("chr_bp" = "chrpos"))
      
      filtered_result
    } else {
      # 如果没有找到匹配的 File_Name，返回空数据框
      data.frame()
    }
  }) %>%
  ungroup()
result_filtered<-as.data.frame(result_filtered)  
write.table(result_filtered,'908-3_ld_0.6_lead_snps_merge_gwas.txt', quote=F ,sep = "\t", row.names = FALSE)

#合并显著的eqtl效应
result_filtered2<-subset(result_filtered,result_filtered$chr_bp %in% eqtl3$chr_bp) #7459
#多对多进行匹配eqtl效应
result_filtered3 <- result_filtered2 %>%
  inner_join(eqtl3, by = c("chr_bp" = "chr_bp"), relationship = "many-to-many")


write.table(result_filtered3,'908-3_ld_0.6_eqtl_RES_unique_snpNgene_sig_withMHC_merge_GWAS_eqtl_res.txt', quote=F ,sep = "\t", row.names = FALSE)

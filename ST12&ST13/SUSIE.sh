######################质控文件待处理的summary数据的文件,生成susie需要的的parquet文件######################

#!/bin/bash
# 设置输入和输出路径
input_path="./step2_160phe_leadsnp"
output_path="./step3_munge"

# 提取所有以mergesnp.txt结尾的文件
mergesnp_files=($(find "$input_path" -name "*mergesnps.txt"))

# 循环处理每个文件
for file_path in "${mergesnp_files[@]}"; do
    # 提取文件名（不包含路径），去掉后缀
    filename=$(basename -- "$file_path")
    output_filename="${filename%_mergesnp.txt}"

    # 构造输出文件路径
    output_file="${output_path}/${output_filename}.parquet"

    # 执行python脚本
    python ./polyfun-master/munge_polyfun_sumstats.py \
      --sumstats "$file_path" \
      --n 35469 \
      --out "$output_file" \
      --min-info 0.3 \
      --min-maf 0.01
done

ls ./step3_munge/*.parquet | wc -l #804，不包括mhc区域



############################################SUSIE############################################

#!/bin/bash
# 设置文件夹路径
input_dir="./step3_munge"
output_dir="./step4_susie/result_susie"

# 循环处理每个以.parquet结尾的文件
for file in "$input_dir"/*.parquet; do
    # 提取文件名（不包含路径）
    filename=$(basename "$file")
    
    # 使用"_"进行拆分
    IFS='_' read -ra parts <<< "$filename"
    
    # 提取所需部分
    type="${parts[0]}"
    type2="${parts[1]}"
    feature="${parts[2]}"
    lobe="${parts[3]}"
    a="${parts[4]}"
    
    # 使用"."进行拆分
    IFS='.' read -ra a_parts <<< "$a"
    
    CHR="${a_parts[0]}"
    START="${a_parts[1]}"
    END2="${a_parts[2]}"
    
    # 运行您的python命令
    python ./polyfun-master/finemapper.py \
        --geno "./step1_extract_160phe_plinkfile/CHR${CHR}.${START}.${END2}" \
        --sumstats "$file" \
        --n 35469 \
        --chr "$CHR" \
        --start "$START" \
        --end "$END2" \
        --method susie \
        --allow-missing \
        --non-funct \
        --max-num-causal 5 \
        --out "${output_dir}/finemap.${type}_${type2}_${feature}_${lobe}_${CHR}.${START}.${END2}.gz"
done














#test,修改py脚本，CS为0.99
input_dir="./step3_munge"
output_dir="./step4_susie/0.99"

# 循环处理每个以.parquet结尾的文件
for file in "$input_dir"/*.parquet; do
    # 提取文件名（不包含路径）
    filename=$(basename "$file")
    
    # 使用"_"进行拆分
    IFS='_' read -ra parts <<< "$filename"
    
    # 提取所需部分
    type="${parts[0]}"
    type2="${parts[1]}"
    feature="${parts[2]}"
    lobe="${parts[3]}"
    a="${parts[4]}"
    
    # 使用"."进行拆分
    IFS='.' read -ra a_parts <<< "$a"
    
    CHR="${a_parts[0]}"
    START="${a_parts[1]}"
    END2="${a_parts[2]}"
    
    # 运行您的python命令
    python /home/jichen/software/polyfun-master/finemapper2.py \
        --geno "./step1_extract_160phe_plinkfile/CHR${CHR}.${START}.${END2}" \
        --sumstats "$file" \
        --n 35469 \
        --chr "$CHR" \
        --start "$START" \
        --end "$END2" \
        --method susie \
        --allow-missing \
        --non-funct \
        --max-num-causal 5 \
        --out "${output_dir}/finemap.${type}_${type2}_${feature}_${lobe}_${CHR}.${START}.${END2}.gz"
done


########################################没跑出来的设置参数，去掉1kg里重复的点，重跑（CS=0.99）##############################
#!/bin/bash

# 指定输入文件和输出目录，input_file只要放没跑出来的parquet文件路径
input_file="./step4_susie/test.txt2" 
output_dir="./step4_susie/0.99"

# 读取输入文件的每一行路径
while IFS= read -r sumstats_file; do
    # 提取文件名（不包含路径）
    filename=$(basename "$sumstats_file")
    
    # 使用"_"进行拆分
    IFS='_' read -ra parts <<< "$filename"
    
    # 提取所需部分
    type="${parts[0]}"
    type2="${parts[1]}"
    feature="${parts[2]}"
    lobe="${parts[3]}"
    a="${parts[4]}"
    
    # 使用"."进行拆分
    IFS='.' read -ra a_parts <<< "$a"
    
    CHR="${a_parts[0]}"
    START="${a_parts[1]}"
    END2="${a_parts[2]}"
    
    # 运行 Python 脚本
    python /home/jichen/software/polyfun-master/finemapper2.py \
        --geno "./step1_extract_160phe_plinkfile/CHR${CHR}.${START}.${END2}" \
        --sumstats "$sumstats_file" \
        --n 35469 \
        --chr "$CHR" \
        --start "$START" \
        --end "$END2" \
        --method susie \
        --allow-missing \
        --non-funct \
        --max-num-causal 5 \
        --susie-resvar 0.9999 \
        --out "${output_dir}/finemap.${type}_${type2}_${feature}_${lobe}_${CHR}.${START}.${END2}.gz"
done < "$input_file"





###############所有位点，其实不用，因为只要credible不是0的才是真正的因果变异！！！！！！！！！！！###############
all=fread('/data1/Imaging_Assoc_2023/calculation_result_JC/SusiE/step4_susie/result_susie/all_traits_result_susie.txt')
library(dplyr)

# 对数据框基于CHR和BP列进行去重，但保留其他列
distinct_all <- distinct(all, CHR, BP, .keep_all = TRUE)
anno_all<-distinct_all[,c(2,3,3,4,5)]
names(anno_all)<-c("CHR","START","END","REF","ALT")

write.table(anno_all, "/data1/Imaging_Assoc_2023/calculation_result_JC/SusiE/step4_susie/susie_annovar/susie_result_160trait_credible_set_snps_for_annovar_ALL.txt", quote=F ,sep = "\t", row.names = FALSE)


cd /data1/Imaging_Assoc_2023/calculation_result_JC/SusiE/step4_susie/susie_annovar
#基本注释
perl /Public/zhuxia/tools/annovar/table_annovar.pl \
    susie_result_160trait_credible_set_snps_for_annovar_ALL.txt \
    /Public/zhuxia/tools/annovar/humandb/ \
    -buildver hg19 \
    -out susie_result_160trait_credible_set_snps_annotion_ALL \
    -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -csvout -polish




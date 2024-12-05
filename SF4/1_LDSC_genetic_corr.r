###########肺叶内部的ldsc#######

#UR.sh 
#same for other 4 lung lobes
#!/bin/bash
cd ./step1_sumstats_result
# 设置路径变量
sumstats_dir="./step1_sumstats_result"
ldscores_dir="./ldsc/ldsc-master/eas_ldscores"
output_dir="./step4_lobe_inner_ldsc/result/UR"

# 获取所有以 _UR_munge.sumstats.gz 结尾的文件列表，并按文件名排序
files=$(ls -1 ${sumstats_dir}/*_UR_munge.sumstats.gz | sort)

# 遍历文件列表
for file1 in $files; do
    for file2 in $files; do
        # 提取文件名（不包括路径和扩展名）
        filename1=$(basename -- "$file1")
        filename2=$(basename -- "$file2")

        # 创建最简单的名字a，b
        a=$(basename -- "$filename1" "_munge.sumstats.gz")
        b=$(basename -- "$filename2" "_munge.sumstats.gz")

        # 确保不重复计算和排除相同文件对
        if [[ "$a" != "$b" && "$file1" < "$file2" ]]; then
            # 执行 ldsc.py 计算遗传相关性
            ./ldsc-master/ldsc.py \
            --rg $filename1,$filename2 \
            --ref-ld-chr ${ldscores_dir}/ \
            --w-ld-chr ${ldscores_dir}/ \
            --out ${output_dir}/${a}.${b}_rg
        fi
    done
done

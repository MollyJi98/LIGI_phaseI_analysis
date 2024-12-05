# LD r01 & LD r08

# r08
for chr in {1..22}
do 
echo ${chr}
plink --bfile ./garfield/chr${chr}_1kgv3_2015_EAS \
--show-tags all \
--tag-kb 500 \
--tag-r2 0.8 \
--out chr${chr}_tag_pri
done

# r01
for chr in {1..22}
do 
echo ${chr}
plink --bfile ./garfield/chr${chr}_1kgv3_2015_EAS \
--show-tags all \
--tag-kb 500 \
--tag-r2 0.1 \
--out chr${chr}_tag_r01_pri
done


# r08
for chr in {1..22}
do 
echo ${chr}
awk '{if($8 != "NONE"){print $3,$8}}' chr${chr}_tag_pri.tags.list > chr${chr}_tag_process_1

sed "s/|/,/g" chr${chr}_tag_process_1 > chr${chr}_tag_process_2
done


R
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
for(i in 1:22){
  re <- fread(paste('chr',i,'_tag_process_2',sep=''),data.table=F)
  re1 <- re[order(re$BP),]
  write.table(re1,file=paste('chr',i,sep=''),row.names=F,col.names=F,quote=F)
}

rm chr*_tag_process_1
rm chr*_tag_process_2


# r01
for chr in {1..22}
do 
echo ${chr}
awk '{if($8 != "NONE"){print $3,$8}}' chr${chr}_tag_r01_pri.tags.list > chr${chr}_tag_process_1

sed "s/|/,/g" chr${chr}_tag_process_1 > chr${chr}_tag_process_2
done

R
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
for(i in 1:22){
  re <- fread(paste('chr',i,'_tag_process_2',sep=''),data.table=F)
  re1 <- re[order(re$BP),]
  write.table(re1,file=paste('chr',i,sep=''),row.names=F,col.names=F,quote=F)
}

rm chr*_tag_process_1
rm chr*_tag_process_2



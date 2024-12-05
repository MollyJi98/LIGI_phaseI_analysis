plink \
--bfile ./chr_all_1kgv3_2015_EAS \
--maf 0.0001 \
--hwe 0.001 \
--geno 0.05 \
--make-bed \
--out chr_all_1kgv3_2015_EAS_MAF0001

# remove duplicated variants
R
options(stringsAsFactors=F)
library(data.table)
setwd('./')
rm(list=ls())
bim <- fread('chr_all_1kgv3_2015_EAS_MAF0001.bim',header = F,data.table=F)
dim(bim)

bim$chr_pos <- paste(bim$V1, bim$V4, sep=':') 
dup <- subset(bim, duplicated(chr_pos)) 
snp <- subset(bim, chr_pos %in% dup$chr_pos) # 28156
nrow(snp)
head(snp)
write.table(snp$V2, file='SNP_rsid_Duplicate.chr_pos',row.names=F,col.names=F,quote=F)

plink \
--bfile chr_all_1kgv3_2015_EAS_MAF0001 \
--exclude SNP_rsid_Duplicate.chr_pos \
--make-bed \
--out chr_all_1kgv3_2015_EAS_mvdup

cp chr_all_1kgv3_2015_EAS_mvdup.bim chr_all_1kgv3_2015_EAS_mvdup.bim.old

# check
R
library(data.table)
rm(list=ls())
bim <- fread('chr_all_1kgv3_2015_EAS_mvdup.bim.old',h=F,data.table=F)
write.table(bim[,c(1,4,3,4,5,6)], file='chr_all_1kgv3_2015_EAS_use.bim',row.names=F,col.names=F,quote=F)

rm(list=ls())
a <- fread('chr_all_1kgv3_2015_EAS_mvdup.bim',h=F,data.table=F)
b <- fread('chr_all_1kgv3_2015_EAS_mvdup.bim.old',h=F,data.table=F)
table(a$V1 == b$V1)
table(a$V2 == b$V2)
table(a$V3 == b$V3)
table(a$V4 == b$V4)
table(a$V5 == b$V5)
table(a$V6 == b$V6)

cp chr_all_1kgv3_2015_EAS_use.bim chr_all_1kgv3_2015_EAS_mvdup.bim

for chr in {1..22}
do 
echo ${chr}

plink \
--bfile chr_all_1kgv3_2015_EAS_mvdup \
--chr ${chr} --make-bed \
--out chr${chr}_1kgv3_2015_EAS
done







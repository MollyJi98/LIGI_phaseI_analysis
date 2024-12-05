# 39服务器
ref_dir='./garfield'

# calculate MAF
for chr in {1..22}
do 
echo ${chr}

plink \
--bfile ${ref_dir}/chr_all_1kgv3_2015_EAS_mvdup \
--chr ${chr} \
--freq \
--out ./chr${chr}
done


# calculate distance to transcription start site (TSS)
R
options(stringsAsFactors=F)
library(data.table)
setwd('./')
rm(list=ls())
bim <- fread('./garfield/chr_all_1kgv3_2015_EAS_mvdup.bim',header = F,data.table=F)
dim(bim)
bim$chr <- paste0('chr',bim$V1)
bim$xuhao <- 1:nrow(bim)
bim$bxuhao <- paste0('b',bim$xuhao)
bim$strand <- '+'
write.table(bim[,c('chr','V4','V4','bxuhao','xuhao','strand','V2')],file='chr_all_1kgv3_2015_EAS_mvdup.bed',row.name=F,col.names=F,quote=F,sep='\t')
system('wc -l chr_all_1kgv3_2015_EAS_mvdup.bed')
system('head chr_all_1kgv3_2015_EAS_mvdup.bed')

# gencode ref
R 
options(stringsAsFactors=F)
library(data.table)
library(stringr)
rm(list=ls())
gtf <- fread("gencode.v19.annotation.gtf", skip = 5)
gtf <- subset(gtf, V3 == "gene")
gtf[, "gene_id" := tstrsplit(V9, ";")[1]]
gtf$gene_id <- gsub('gene_id "', '', gtf$gene_id)
gtf$gene_id <- gsub('"', '', gtf$gene_id)
gtf[, "gene_name" := tstrsplit(V9, ";")[5]]
gtf$gene_name <- gsub(' gene_name "', '', gtf$gene_name)
gtf$gene_name <- gsub('"', '', gtf$gene_name)
gtf <- gtf[,c(1,4,5,7,10,11)]
names(gtf)[1:4] <- c('chr','start','end','direction')
write.csv(gtf,file = 'gene.csv',row.names = F,quote = F)

gen <- read.csv('gene.csv')
gen$tss <- ifelse(gen$direction=='+',gen$start,ifelse(gen$direction=='-',gen$end, NA))
gen$chr1 <- gen$chr
gen$chr1 <- as.numeric(gsub('chr','',gen$chr1))
gen <- subset(gen, !is.na(chr1))
gen <- gen[order(gen$chr1,gen$tss),]

gen$xuhao <- 1:nrow(gen)
gen$bxuhao <- paste('a',gen$xuhao,sep='')

write.table(gen[,c('chr','tss','tss','bxuhao','xuhao','direction')],file='garfield.gencode.v19.gene.ann.bed',row.name=F,col.names=F,quote=F,sep='\t')
system('wc -l garfield.gencode.v19.gene.ann.bed')
system('head garfield.gencode.v19.gene.ann.bed')


bedtools closest -a chr_all_1kgv3_2015_EAS_mvdup.bed \
-b garfield.gencode.v19.gene.ann.bed -D b \
> chr_all_1kgv3_2015_EAS_mvdup_TSS_Distance

R
setwd('./garfield/maftssd')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())

# Distance to nearest TSS
dat <- fread('chr_all_1kgv3_2015_EAS_mvdup_TSS_Distance',h=F,data.table=F) 
head(dat)
dat <- subset(dat, !duplicated(V4)) 
dim(dat)
dat$V1 <- gsub('chr','',dat$V1)

# Combine MAF with "Distance to nearest TSS"
for(i in 1:22){
  maf <- fread(paste0('chr',i,'.frq'))
  tss <- subset(dat, gsub('chr','',dat$V1) == i)
  print(paste0('chr',i,' MAF:',nrow(maf),'; TSS:',nrow(tss)))
  print(nrow(maf)==nrow(tss))
  
  re <- merge(maf[,c('SNP','MAF')], tss[,c('V2','V7','V14')], by.x='SNP', by.y='V7')
  re_final <- re[order(re$V2),]
  re_final <- re_final[,c(3,2,4)]
  write.table(re_final, file=paste0('chr',i), row.names=F, col.names=F, quote=F)
  
  print(nrow(re_final)==nrow(tss))
  rm(list=c('maf','tss','re','re_final'))
}



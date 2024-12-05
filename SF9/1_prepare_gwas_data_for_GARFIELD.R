# firstorder
R
options(stringsAsFactors=F)
library(ggplot2)
library(data.table)
rm(list=ls())

data<-fread('firstorder_minP.txt') 
dim(data)  

for(i in 1:22){
  sda <- subset(data, CHR==i)
  sda1 <- sda[order(sda$BP),c('BP','P_BOLT_LMM')]  
  write.table(sda1,file=paste0('./garfield/pval/chr',i),row.names=F,col.names=F,quote=F)
  print(i)
}


# shape
R
options(stringsAsFactors=F)
library(ggplot2)
library(data.table)
rm(list=ls())

data<-fread('shape_minP.txt') 
dim(data)  

for(i in 1:22){
  sda <- subset(data, CHR==i)
  sda1 <- sda[order(sda$BP),c('BP','P_BOLT_LMM')]  
  write.table(sda1,file=paste0('./garfield/pval/chr',i),row.names=F,col.names=F,quote=F)
  print(i)
}

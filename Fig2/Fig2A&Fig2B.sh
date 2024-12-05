#################################################mplot#################################################
R
library(data.table)

common_chr_path <- "./Assoc_Lung_phenotype_low_3.8w_clean"
sex_chr_path <- "./Assoc_Lung_phenotype_low_3.8w_clean_chrX"
common_chr_files <- list.files(common_chr_path, full.names = TRUE)
sex_chr_files <- list.files(sex_chr_path, full.names = TRUE)

output_path <- "./combine_chrX_0.05"

# loop for data extraction
for (i in 1:length(common_chr_files)) {
  merged_data <- data.frame()  
  # 读取常染色体文件并保留 P_BOLT_LMM < 0.05 的行，选择 CHR、BP 和 P_BOLT_LMM 三列
  common_chr_data <- fread(common_chr_files[i])
  common_chr_data <- common_chr_data[P_BOLT_LMM < 0.05, .(CHR, BP, P_BOLT_LMM)]
  names(common_chr_data)<-c('CHR','BP','P')
  sex_chr_data <- fread(sex_chr_files[i])
  sex_chr_data<-subset(sex_chr_data,sex_chr_data$P < 0.05)
  sex_chr_data <- sex_chr_data[,c("#CHROM","POS","P")]
  names(sex_chr_data)<-c('CHR','BP','P')
  
  # merge
  merged_data <- rbind(common_chr_data, sex_chr_data)

  # output
  output_file <- paste0(basename(common_chr_files[i]), "_combine_chrX_0.05.txt")
  write.table(merged_data, file.path(output_path, output_file), sep = "\t", row.names = FALSE,quote=F)
}


############data prepare############
R
library(data.table)
library(ggplot2)

file= list.files(pattern = '_combine_chrX_0.05.txt')     

file2 <- grep('_combine_chrX_0.05.txt', file, value = TRUE)
names <- file2
n = length(file2) #160

dir = paste("./",file2,sep="")  
color<-c(rep(c("#6495ED","#384075"),11),"#6495ED")

remove=fread('./significant_remove_outlier_chr_bp.txt') #outliersSNPs


###first order###
matching_names <- names[grep("^original_firstorder", names)]
all_data <- rbindlist(lapply(matching_names, fread))
names(all_data) <- c("CHR", "BP", "P_BOLT_LMM")
all_data$CHR_BP <- paste(all_data$CHR, all_data$BP, sep = ":")
result <- all_data[order(P_BOLT_LMM), .SD[1], by = CHR_BP]
result <- result[!CHR_BP %in% remove$CHR_BP]
result$CHR[result$CHR == "X"] <- 23

# save
write.table(result,"./firstorder_0.05_snp_for_mplot.txt",quote=F, sep = "\t", row.names = FALSE)


########################plot-firstorder########################
library(data.table)
library(ggplot2)
library(dplyr)

result<-fread('./firstorder_0.05_snp_for_mplot.txt') 
sig_first<-subset(result,result$P_BOLT_LMM<5e-8) #highlight
highlight<-fread("./firstorder_highlight_range.txt") #highlight250kb

#任意选一个表型，挑出实际要highlight的点（highlight上下游250kb）
a=fread('./original_firstorder_10Percentile_CR_assoc.result')
b=fread('./original_firstorder_10Percentile_CR_assoc.lung.low.imputed.chrX.original_firstorder_10Percentile_CR.glm.linear')

b$CHR<-"23"
names(b)[3]<-"BP"
b$chrpos<-paste(b$CHR,b$BP,sep=":")


highlighted_data <- rbindlist(lapply(1:nrow(highlight), function(i) {
  a_subset <- a[CHR == highlight$CHR[i] & BP >= highlight$START[i]-250000 & BP <= highlight$END[i]+250000, "chrpos"]
  return(a_subset)
}))

highlighted_data2 <- rbindlist(lapply(1:nrow(highlight), function(i) {
  b_subset <- b[CHR == highlight$CHR[i] & BP >= highlight$START[i]-250000 & BP <= highlight$END[i]+250000, "chrpos"]
  return(b_subset)
}))

chrpos_vector <- highlighted_data2$chrpos
chrpos_vector <- gsub("23:", "X:", chrpos_vector)
highlighted_data2$chrpos <- chrpos_vector
highlighted_data_final<-rbind(highlighted_data,highlighted_data2) #最终常染色体和性染色体要highlight的内容（显著位点上下游250kb）


dat <- result[!is.na(result$P_BOLT_LMM),]
d=transform(dat,logp=-log10(P_BOLT_LMM))
d$pos<-as.numeric(d$BP)

d2<-subset(d,d$P_BOLT_LMM<0.001) 
df.tmp<- d2 %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(pos)) %>%
  mutate(tot= cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(d2, ., by=("CHR"="CHR")) %>%
  arrange(CHR,pos) %>%
  mutate(BPcum=pos+tot)

# 构建is_highlight标签
df.tmp$is_highlight <- "no"
df.tmp$is_highlight<-ifelse(df.tmp$CHR_BP %in% highlighted_data_final$chrpos,"yes",df.tmp$is_highlight)


# 检查结果
head(df.tmp)  
axisdf<-  df.tmp%>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)

subtmp=subset(df.tmp,is_highlight=="yes")

mypalette<-c("gray70","gray90")
col<-mypalette

pdf("./manhattan_first_250000_0.001.pdf", width = 15, height = 8) 

ggplot(df.tmp, aes(x=BPcum, y=logp)) +
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(col, 23 )[1:23]) +
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks=axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) + # expand-c(0.0)removesspace between plot area and axis
  # add plot and axis titles
  # ggtitle(paste0(title)) +
  labs(x ="Chromosome") +
  labs(y ="-log10(p)")+
  # add genome wide sig and sugg lines
  geom_hline(yintercept = -log10(5e-8), colour="red", linetype = "dashed")+
  # geom_hline(yintercept = -log10(sugg), linetype-"dashed") +
  # Add highlighted points
  geom_point(data= subset(df.tmp, is_highlight=="yes"), color="royalblue", size=0.5) +
  theme_bw(base_size =15) +
  theme(
    plot.title = element_text(hjust = 0.5, size =15),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    axis.line = element_line(color = "black", size = 0.5),
  )


dev.off()


###shape###
matching_names <- names[grep("^original_shape", names)]
all_data <- rbindlist(lapply(matching_names, fread))
names(all_data) <- c("CHR", "BP", "P_BOLT_LMM")
all_data$CHR_BP <- paste(all_data$CHR, all_data$BP, sep = ":")
result <- all_data[order(P_BOLT_LMM), .SD[1], by = CHR_BP]
result <- result[!CHR_BP %in% remove$CHR_BP]

result$CHR[result$CHR == "X"] <- 23

write.table(result,"./shape_0.05_snp_for_mplot.txt",quote=F, sep = "\t", row.names = FALSE)

########################plot-shape########################
library(data.table)
library(ggplot2)
library(dplyr)

result<-fread('./shape_0.05_snp_for_mplot.txt')
sig_shape<-subset(result,result$P_BOLT_LMM<5e-8)
highlight<-fread("./shape_highlight_range.txt")
a=fread('./original_firstorder_10Percentile_CR_assoc.result')
b=fread('./original_firstorder_10Percentile_CR_assoc.lung.low.imputed.chrX.original_firstorder_10Percentile_CR.glm.linear')

b$CHR<-"23"
names(b)[3]<-"BP"
b$chrpos<-paste(b$CHR,b$BP,sep=":")


highlighted_data <- rbindlist(lapply(1:nrow(highlight), function(i) {
  a_subset <- a[CHR == highlight$CHR[i] & BP >= highlight$START[i]-250000 & BP <= highlight$END[i]+250000, "chrpos"]
  return(a_subset)
}))

highlighted_data2 <- rbindlist(lapply(1:nrow(highlight), function(i) {
  b_subset <- b[CHR == highlight$CHR[i] & BP >= highlight$START[i]-250000 & BP <= highlight$END[i]+250000, "chrpos"]
  return(b_subset)
}))


chrpos_vector <- highlighted_data2$chrpos
chrpos_vector <- gsub("23:", "X:", chrpos_vector)
highlighted_data2$chrpos <- chrpos_vector

highlighted_data_final<-rbind(highlighted_data,highlighted_data2)

dat <- result[!is.na(result$P_BOLT_LMM),]

d=transform(dat,logp=-log10(P_BOLT_LMM))
d$pos<-as.numeric(d$BP)

d2<-subset(d,d$P_BOLT_LMM<0.001)
df.tmp<- d2 %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(pos)) %>%
  mutate(tot= cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(d2, ., by=("CHR"="CHR")) %>%
  arrange(CHR,pos) %>%
  mutate(BPcum=pos+tot)


df.tmp$is_highlight <- "no"
df.tmp$is_highlight<-ifelse(df.tmp$CHR_BP %in% highlighted_data_final$chrpos,"yes",df.tmp$is_highlight)


axisdf<-  df.tmp%>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)
subtmp=subset(df.tmp,is_highlight=="yes")
mypalette<-c("gray70","gray90")
col<-mypalette

pdf("./manhattan_shape_250000_0.001.pdf", width = 15, height = 8) 

ggplot(df.tmp, aes(x=BPcum, y=logp)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(col, 23 )[1:23]) +
  scale_x_continuous(label = axisdf$CHR, breaks=axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  labs(x ="Chromosome") +
  labs(y ="-log10(p)")+
  geom_hline(yintercept = -log10(5e-8), colour="red", linetype = "dashed")+
  geom_point(data= subset(df.tmp, is_highlight=="yes"), color="royalblue", size=0.5) +
  theme_bw(base_size =15) +
  theme(
    plot.title = element_text(hjust = 0.5, size =15),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    axis.line = element_line(color = "black", size = 0.5)
  )


dev.off()

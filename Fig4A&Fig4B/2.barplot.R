#################s-ldsc plot#################
library(tidyverse)
library(RColorBrewer)
library(dendextend)
library(tibble)
library(ggplotify)
library(data.table)
library(dplyr)
library(ggrepel)
library(Cairo)
library(grid)

# shape
ldsc.dat <- fread("shape_minP_s_ldsc_res.txt")
ldsc.dat <-subset(ldsc.dat,ldsc.dat$Plot==1)
ldsc.dat$names <- paste0(ldsc.dat$`Chromatin mark`,' (',ldsc.dat$Category,')')
ldsc.dat$p.log10 <- -log10(ldsc.dat$`Enrichment P-value`)
ldsc.final <- ldsc.dat
#ldsc.final$cells.ordered<- factor(ldsc.final$`Chromatin mark`, c("DNase","H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac"))

small_ldsc_theme <- theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust = 1.1,
                                   size = 8,family = "Arial",color="black",face="bold"),
        legend.position ="right",
        plot.title = element_text(hjust = 0.5, face="bold",size=10,family="Arial", vjust = -1),
        axis.text.y = element_text(size = 6,colour = "black",family = "Arial",face="bold"),
        axis.title.x = element_text(size=6,color = "black",family = "Arial",face = "bold"),
        axis.title.y = element_text(size=8,color = "black",family = "Arial",face = "bold"),
        #panel.grid.major.y = element_line(color = "grey80",linetype = "dashed"),
        #plot.title = element_blank(),
        legend.text = element_text(size = 6,family = "Arial",color = "black"),
        legend.spacing.y = unit(0.2,'cm'),
        legend.spacing.x = unit(0.2,'cm'),
        legend.background = element_blank(),
        legend.key.size = unit(1, 'mm'))

col <- c(`Aorta`="#4DAF4A",
         `Brain Anterior Caudate`="#377EB8",
         `Fetal Lung`="#984EA3",
         `Lung`= "#E41A1C",
         `Esophagus`="#2D417F",
         `Gastric`="#FF7F00",
         `Left Ventricle`="#F781BF",
         `Liver`="#8c564b",
         `Pancreas`="#7f7f7f",
         `Rectal Mucosa Donor 29`="#bcbd22",
         `Small Intestine`="#17becf",
         `Spleen`="#aec7e8",
         `Thymus`="#ffbb78",
         `Colonic Mucosa`="#98df8a"
)

ldsc.final$type <- factor(ldsc.final$names)
ldsc.final$plot.order <-  factor(ldsc.final$Category,levels = c('Fetal Lung','Lung','Aorta','Brain Anterior Caudate',
                                                                'Esophagus','Gastric','Left Ventricle','Liver',
                                                                'Pancreas','Rectal Mucosa Donor 29','Small Intestine',
                                                                'Spleen','Thymus','Colonic Mucosa'))

p.min <- 0
p.max <- 25
p.sig <- -log10(0.05/(145))

ldsc.final2 <- ldsc.final %>%
  arrange(plot.order)

ldsc.final2 <- ldsc.final2[-43,] #remove invalid value


##########拼图##############
ldsc.final2$plot.order <- factor(ldsc.final2$plot.order, levels = unique(ldsc.final2$plot.order))  

# 过滤掉 p.log10 为 NA 的行  
ldsc.final2 <- ldsc.final2 %>%  
  filter(!is.na(p.log10)) %>%  
  arrange(`Chromatin mark`, desc(p.log10))  

# 使用 facet_grid 垂直分组显示每个 Chromatin mark  
p.bar <- ggplot(ldsc.final2, aes(  
  x = p.log10,                
  y = reorder(type, p.log10),  
  fill = Category)) +        
  geom_bar(stat = "identity", width = 0.8) +  
  geom_vline(xintercept = p.sig, color = "blue", linetype = "longdash", size = 0.75, alpha = 0.75) + 
  
  # 设置标题和轴标签  
  ggtitle("Shape") +   
  ylab("Category") + 
  xlab("-log10(Coefficient P-value)") +  
  
  scale_x_continuous(limits = c(0, p.max), expand = c(0.01, 0)) +   
  scale_fill_manual(values = col) +   
  
  # 添加主题
  small_ldsc_theme +   
  theme(legend.position = "right",   
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  
        axis.text.x = element_text(angle = 0, hjust = 0.4),  # 旋转 x 轴标签以便阅读  
        axis.text.y = element_text(margin = margin(0, 0, 0, 0, "cm"))) +   
  # 使用 facet_grid 垂直分组显示每个 Chromatin mark  
  facet_grid(rows = vars(`Chromatin mark`), scales = "free_y", space = "free_y")+   
  theme(
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#f3f3f3", colour = "NA"),
    strip.text = element_text(size = 8,family = "Arial",color="black",face="bold"))

print(p.bar)
# 显示图表  
CairoPDF("shape_s_ldsc_barplot.pdf", height = 6, width = 5)
print(p.bar)
dev.off()



# firstorder
ldsc.dat <- fread("firstorder_minP_s_ldsc_res.txt")
ldsc.dat <-subset(ldsc.dat,ldsc.dat$Plot==1)
ldsc.dat$names <- paste0(ldsc.dat$`Chromatin mark`,' (',ldsc.dat$Category,')')
ldsc.dat$p.log10 <- -log10(ldsc.dat$`Enrichment P-value`)
ldsc.final <- ldsc.dat


small_ldsc_theme <- theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust = 1.1,
                                   size = 8,family = "Arial",color="black",face="bold"),
        legend.position ="right",
        plot.title = element_text(hjust = 0.5, face="bold",size=10,family="Arial", vjust = -1),
        axis.text.y = element_text(size = 6,colour = "black",family = "Arial",face="bold"),
        axis.title.x = element_text(size=6,color = "black",family = "Arial",face = "bold"),
        axis.title.y = element_text(size=8,color = "black",family = "Arial",face = "bold"),
        legend.text = element_text(size = 6,family = "Arial",color = "black"),
        legend.spacing.y = unit(0.1,'cm'),
        legend.spacing.x = unit(0.1,'cm'),
        legend.background = element_blank(),
        legend.key.size = unit(1, 'mm'))

col <- c(`Aorta`="#4DAF4A",
         `Brain Anterior Caudate`="#377EB8",
         `Fetal Lung`="#984EA3",
         `Lung`= "#E41A1C",
         `Esophagus`="#2D417F",
         `Gastric`="#FF7F00",
         `Left Ventricle`="#F781BF",
         `Liver`="#8c564b",
         `Pancreas`="#7f7f7f",
         `Rectal Mucosa Donor 29`="#bcbd22",
         `Small Intestine`="#17becf",
         `Spleen`="#aec7e8",
         `Thymus`="#ffbb78",
         `Colonic Mucosa`="#98df8a"
)

ldsc.final$type <- factor(ldsc.final$names)
ldsc.final$plot.order <-  factor(ldsc.final$Category,levels = c('Fetal Lung','Lung','Aorta','Brain Anterior Caudate',
                                                                'Esophagus','Gastric','Left Ventricle','Liver',
                                                                'Pancreas','Rectal Mucosa Donor 29','Small Intestine',
                                                                'Spleen','Thymus','Colonic Mucosa'))

p.min <- 0
p.max <- 10
p.sig <- -log10(0.05/(145))

ldsc.final2 <- ldsc.final %>%
  arrange(plot.order)


##########多个柱状图拼图##############
ldsc.final2$plot.order <- factor(ldsc.final2$plot.order, levels = unique(ldsc.final2$plot.order))  

# 过滤掉 p.log10 为 NA 的行  
ldsc.final2 <- ldsc.final2 %>%  
  filter(!is.na(p.log10)) %>%  
  arrange(`Chromatin mark`, desc(p.log10))  

# 重新生成柱状图，使用 facet_grid 垂直分组显示每个 Chromatin mark  
p.bar <- ggplot(ldsc.final2, aes(  
  x = p.log10,                
  y = reorder(type, p.log10), 
  fill = Category)) +             
  
  geom_bar(stat = "identity", width = 0.8) +  
  geom_vline(xintercept = p.sig, color = "blue", linetype = "longdash", size = 0.75, alpha = 0.75) + 
  
  ggtitle("First order") +   
  ylab("Category") +  
  xlab("-log10(Coefficient P-value)") +   
  

  scale_x_continuous(limits = c(0, p.max), expand = c(0.01, 0)) +   
  scale_fill_manual(values = col) +   
  
  small_ldsc_theme +   
  theme(legend.position = "right",   
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  
        axis.text.x = element_text(angle = 0, hjust = 0.4),  # 旋转 x 轴标签以便阅读  
        axis.text.y = element_text(margin = margin(0, 0, 0, 0, "cm"))) +   
  # 使用 facet_grid 垂直分组显示每个 Chromatin mark  
  facet_grid(rows = vars(`Chromatin mark`), scales = "free_y", space = "free_y") +  
  theme(
    legend.title = element_blank(),
   strip.background = element_rect(fill = "#f3f3f3", colour = "NA"),
   strip.text = element_text(size = 8,family = "Arial",color="black",face="bold"))

# 显示图表  
print(p.bar)

CairoPDF("firstorder_s_ldsc_barplot.pdf", height = 6, width = 5)
print(p.bar)
dev.off()



#######自己拼图美化#######




rm(list=ls())

# BiocManager::install("org.Hs.eg.db") #人全基因组注释R包
# BiocManager::install("clusterProfiler")
# install.packages("ggplot2")

#加载所需要的包
library("clusterProfiler")
library(org.Hs.eg.db)
library("ggplot2")
library(openxlsx)


#############################加gwama结果##########################
#文件预处理
inputFile=read.xlsx("两类_near_gene_coloc整理.xlsx",sheet = "GWAMA_sig")
first <-  inputFile$first_all
shape <-  inputFile$shape_all

###################################################################################
#######################################first#######################################
###################################################################################
#将Symbol转换为EntrezID
diff_entrez<-bitr(
  first,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Hs.eg.db'
)
head(diff_entrez)


#######################################first-GO################################################
go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
#将RNTREZ转换为Symbol
go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')

#提取goG富集结果表格 
go_result<-go_enrich@result
go_result

write.table(go_result,"GO_RES_first_final_with_GWAMA3.xlsx", quote=F ,sep = "\t", row.names = FALSE)


###################################################################################
#######################################shape#######################################
###################################################################################
#将Symbol转换为EntrezID
diff_entrez<-bitr(
  shape,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Hs.eg.db'
)
head(diff_entrez)


#######################################shape-GO################################################
go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
#将RNTREZ转换为Symbol
go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')

#提取goG富集结果表格 
go_result<-go_enrich@result
go_result

write.table(go_result,"GO_RES_shape_final_with_GWAMA3.xlsx", quote=F ,sep = "\t", row.names = FALSE)





############################绘制气泡图#####################
library(ggplot2)
library(openxlsx)
first_go <- read.xlsx('/GO_RES_first_final_with_GWAMA3.xlsx')

# 按照p值大小排序
df <- first_go[order(first_go$p.adjust), ]
df <- subset(df,df$ONTOLOGY=="BP")
# 创建气泡图（纵向）
p2 <- ggplot(df[1:15,], aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust), size = -log10(p.adjust))) +
  geom_point(aes(color = p.adjust), alpha = 0.9) +
  scale_size_continuous(range = c(3, 9)) +
  scale_color_gradient(low = "#5b98e0", high = "firebrick3") +
  labs(title = "GO Pathway Enrichment Analysis of first order traits associated gene",
       x = "-log10(p-adjust)",
       y = "Pathway",
       size = "-log10(p-adjust)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, hjust = 1)) 

p2
ggsave("go_res_plot_firstorder.pdf", plot = p2, device = "pdf",width = 12,height = 4)



shape_go <- read.xlsx('GO_RES_shape_final_with_GWAMA3.xlsx')

# 按照p值大小排序
df <- shape_go[order(shape_go$p.adjust), ]

# 创建气泡图（纵向）
p2 <- ggplot(df[1:15,], aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust), size = -log10(p.adjust))) +
  geom_point(aes(color = p.adjust), alpha = 0.9) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "#5b98e0", high = "firebrick3") +
  labs(title = "GO Pathway Enrichment Analysis of shape traits associated gene",
       x = "-log10(p-adjust)",
       y = "Pathway",
       size = "-log10(p-adjust)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, hjust = 1))
  #coord_cartesian(xlim = c(4.4, 8.6))  # 设置x轴范围
p2
# 保存p2图为PDF文件，指定保存路径
ggsave("go_res_plot_shape.pdf", plot = p2, device = "pdf",width = 12,height = 4)


########venn########
library(openxlsx)
library(ggvenn)
library(eulerr)

####################################成比例文氏图############################################
# 读取并整理  数据文件
venn_dat  <-  read.xlsx("venn_plot_data.xlsx") #需要画图的基因

venn_list <- as.list(venn_dat)                         # 转换成列表
venn_list <- purrr::map(venn_list, na.omit)            # 删除列表中每个向量中的NA
venn_list <- lapply(venn_list, function(x) x[x != ""]) # 删除列表中每个向量中的""空字符串
venn_list <- lapply(venn_list, unique)                 # 移除重复元素
venn_list

pdf("venn_plot.pdf", width = 8, height = 6) 

# plot
plot(euler(
  venn_list,
  shape = "circle"),                    # 图案的形状，椭圆ellipse 或圆circle
  quantities = list(type = c("percent","counts"),cex=1),          # 显示类型，百分比和数字，数字大小
  labels=list(cex=1),                   # 组名标签的大小
  edges = list(col = "black", lex = 2), # 图形边缘的颜色和大小
  fills = list(fill = c("#f7c9c1","#c9eaf2","#f9f8da"),alpha=0.7) # 填充的颜色和透明度
  # legend = list(side = "right")       # 图例的位置
)

dev.off()

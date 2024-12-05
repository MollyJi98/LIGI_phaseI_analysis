#############画四种基因mapping方式的韦恩图#############
#first order和shape
library(openxlsx)
library(ggvenn)
weien <- read.xlsx("两类_near_gene_coloc整理.xlsx",sheet="GWAMA_sig")
#提取firstorder和shape的特征
first <- weien[,c(1:4)]
names(first) <- c('MAGMA','posmap','coloc','SuSIE')
shape <- weien[,c(6:9)]
names(shape) <- c('MAGMA','posmap','coloc','SuSIE')

#first#
# 读取数据文件
venn_dat  <- first
venn_list <- as.list(venn_dat)              # 制作韦恩图搜所需要的列表文件
venn_list <- purrr::map(venn_list, na.omit) # 删除列表中每个向量中的NA

# 绘图
ggvenn(
  data = venn_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#f7c9c1", "#c9eaf2", "#b2e2db", "#c5cbdb"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)

# 保存图片
ggsave("韦恩图_firstorder.pdf", width = 8, height = 6)



#shape#
# 读取数据文件
venn_dat  <- shape
venn_list <- as.list(venn_dat)              # 制作韦恩图搜所需要的列表文件
venn_list <- purrr::map(venn_list, na.omit) # 删除列表中每个向量中的NA

# 绘图
ggvenn(
  data = venn_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#f7c9c1", "#c9eaf2", "#b2e2db", "#c5cbdb"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)

# 保存图片
ggsave("韦恩图_shape.pdf", width = 8, height = 6)

#画完之后ai拼图
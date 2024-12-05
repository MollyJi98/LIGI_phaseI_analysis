#################画位点的upset图#################
#devtools::install_github("krassowski/complex-upset")
rm(list=ls())
library(openxlsx)
library(ComplexUpset)
library(ggplot2)

data_first <- read.xlsx("first_shape_7表型显著情况_for_upset.xlsx",sheet = 1)
upset_first <- data_first[,c(10:16)]
data_shape <- read.xlsx("first_shape_7表型显著情况_for_upset.xlsx",sheet = 2)
upset_shape <- data_shape[,c(10:16)]

###################整合两类表型########################
upset_first$type <- "firstorder"
upset_shape$type <- "shape"
data_combine <- rbind(upset_first,upset_shape)

upset_data <- data_combine[rowSums(data_combine[,1:7]) > 0, ]

# 生成 upset 图，按 type 不同颜色区分
upset_plot <- upset(
  upset_data,
  intersect = colnames(upset_data)[1:7],  # 使用前三列表示交集
  width_ratio = 0.7,
  sort_intersections_by = "cardinality", 
  stripes = upset_stripes(color = c("#f6f6f6", "white")),
  base_annotations = list(
    'Intersection size' = (
      intersection_size(
        counts = TRUE,
        aes(fill = type)  # 按 type 映射颜色
      ) + 
          scale_fill_manual(values = c("firstorder" = "firebrick3", "shape" = "#5b98e0")) +  # 定义颜色
        theme_minimal()
    )
  )
)  + 
  # 修改点的颜色和大小
  scale_color_manual(values = c("gray","#000c7d")) + # 自定义颜色
  scale_size_manual(values = c(1,1)) + # 设置点大小
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # 去除x轴文字
        axis.ticks.x = element_blank())  # 去除x轴刻度线

# 展示图形
upset_plot

ggsave(upset_plot, filename = "ComplexUpset_combine.pdf",
       width = 15, height = 4, device = cairo_pdf)


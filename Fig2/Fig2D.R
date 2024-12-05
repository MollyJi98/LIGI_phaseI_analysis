#first order and shape comparison#
library("ggplot2")
library("ggsci")
library(showtext)
library(Cairo)
library(grid)


data=read.csv("first_order_shape_Compare_v1.csv") #收集整理待比较区域的leadsnp的效应

legend <- ggplot(data) + 
  geom_point(aes(x = BETA_shape, y = BETA, color = Gene)) + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  # 设置透明背景
        plot.background = element_rect(fill = "transparent", colour = NA),  # 设置透明背景
        axis.line = element_line(colour = "#323232", size = 0.4),  # 坐标轴线颜色
        axis.ticks = element_line(color = "#323232", size = 0.2)) +  # 坐标轴刻度线颜色
  scale_x_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.1)) + 
  scale_y_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.1)) + 
  geom_segment(aes(y = BETA_IL, yend = BETA_UL, x = BETA_shape), size = 0.3, color = "#323232") +  # 设置线条颜色
  geom_segment(aes(y = BETA, x = BETA_shape_IL, xend = BETA_shape_UL), size = 0.3, color = "#323232") +  # 设置线条颜色
  geom_vline(xintercept = 0, size = 0.4, color = "#323232") +  # 竖线颜色
  geom_hline(yintercept = 0, size = 0.4, color = "#323232") +  # 横线颜色
  labs(x = "Per-allele effect size for shape features", 
       y = "Per-allele effect size for first order features", 
       color = "Gene") + 
  scale_color_rickandmorty() + 
  theme(panel.grid = element_blank(),  # 去掉背景网格线
        axis.title = element_text(face = "plain", color = "#323232", family = "Arial"),  # 坐标轴标题颜色
        axis.text = element_text(face = "plain", color = "#323232", family = "Arial"),  # 坐标轴文本颜色
        plot.caption = element_text(face = "plain", color = "#323232", family = "Arial"))  # 图注颜色

legend

CairoPDF("first_order_shape_Compare_legend.pdf", height = 4, width = 5)
print(legend)
dev.off()

########################剩余在ppt拼图########################
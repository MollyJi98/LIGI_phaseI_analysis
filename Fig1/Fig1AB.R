R
rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(grid) 
library(svglite)

firstorder <- read_excel("heri.xlsx",sheet="firstorder")
shape <- read_excel("heri.xlsx",sheet="shape")

firstorder <- firstorder %>%
  mutate(Lobe = factor(Lobe, levels = c("UR", "LR", "CR", "UL", "LL"))) %>% 
  arrange(var, Lobe)

shape <- shape %>%
  mutate(Lobe = factor(Lobe, levels = c("UR", "LR", "CR", "UL", "LL"))) %>% 
  arrange(var, Lobe)

### ### ### ### ### ### ### ###firstorder
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
order_firstorder <- c("Energy", "TotalEnergy", "Maximum", "Range", "Kurtosis", "Skewness", 
                      "Uniformity", "Minimum", "10Percentile", "90Percentile", "Entropy", 
                      "InterquartileRange", "MeanAbsoluteDeviation", "Mean", "Median", 
                      "RobustMeanAbsoluteDeviation", "RootMeanSquared", "Variance") 
p_firstorder <- ggplot(firstorder, aes(x = factor(var, levels = order_firstorder), y = her)) +
  geom_violin(trim = FALSE, fill = "grey90", color = "grey90") +
  geom_jitter(aes(color = Lobe, shape = Lobe), width = 0.2, size = 1) + 
  scale_color_manual(values = c("UR" = adjustcolor("#377EB8", alpha.f = 0.7), 
                                "UL" = adjustcolor("#E41A1C", alpha.f = 0.7), 
                                "LR" = adjustcolor("#377EB8", alpha.f = 0.7), 
                                "LL" = adjustcolor("#E41A1C", alpha.f = 0.7), 
                                "CR" = adjustcolor("#377EB8", alpha.f = 0.7))) +  
  scale_shape_manual(values = c("UL" = 16, "UR" = 16, "CR" = 17, "LL" = 15, "LR" = 15)) +  
  labs(x = "Phenotype", y = "Heritability") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(linewidth = 0.15, color = "black"),
    axis.line.x = element_line(arrow = arrow(type = "closed", length = unit(0.15, "cm"))),
    axis.line.y = element_line(arrow = arrow(type = "closed", length = unit(0.15, "cm"))),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 6, color = "black"), 
    axis.title.x = element_text(size = 7, color = "black"),  
    axis.title.y = element_text(size = 7, color = "black"),  
    legend.position = "bottom" 
  )
ggsave(filename = "firstorder.pdf", plot = p_firstorder, width = 5, height = 4)


### ### ### ### ### ### ### ###shape
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
order_shape <- c( "SurfaceVolumeRatio", "Elongation", "Flatness", "Sphericity", 
                  "LeastAxisLength", "MajorAxisLength", "Maximum2DDiameterColumn", 
                  "Maximum2DDiameterRow", "Maximum2DDiameterSlice", "Maximum3DDiameter", 
                  "MeshVolume", "MinorAxisLength", "SurfaceArea", "VoxelVolume") 
p_shape <- ggplot(shape, aes(x = factor(var, levels = order_shape), y = her)) +
  geom_violin(trim = FALSE, fill = "grey90", color = "grey90") +
  geom_jitter(aes(color = Lobe, shape = Lobe), width = 0.2, size = 1) +  
  scale_color_manual(values = c("UR" = adjustcolor("#377EB8", alpha.f = 0.7), 
                                "UL" = adjustcolor("#E41A1C", alpha.f = 0.7), 
                                "LR" = adjustcolor("#377EB8", alpha.f = 0.7), 
                                "LL" = adjustcolor("#E41A1C", alpha.f = 0.7), 
                                "CR" = adjustcolor("#377EB8", alpha.f = 0.7))) +  
  scale_shape_manual(values = c("UL" = 16, "UR" = 16, "CR" = 17, "LL" = 15, "LR" = 15)) +  
  labs(x = "Phenotype", y = "Heritability") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(linewidth = 0.15, color = "black"),  
    axis.line.x = element_line(arrow = arrow(type = "closed", length = unit(0.15, "cm"))),
    axis.line.y = element_line(arrow = arrow(type = "closed", length = unit(0.15, "cm"))),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "black"), 
    axis.text.y = element_text(size = 6, color = "black"),  
    axis.title.x = element_text(size = 7, color = "black"), 
    axis.title.y = element_text(size = 7, color = "black"),   
    legend.position = "bottom"  
  )

ggsave(filename = "shape.pdf", plot = p_shape, width = 5, height = 4)







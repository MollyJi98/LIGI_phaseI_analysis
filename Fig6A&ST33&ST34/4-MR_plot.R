##############JC改配色，源代码来自myl##############
rm(list = ls())

library(tidyverse)
library(data.table)
library(openxlsx)
library(rio)
# 需要确保 0.9.5版本的 ggrepel，不然就要下载后重装!!!
# remove.packages("ggrepel" )
# install.packages('/Users/molly/Downloads/ggrepel-0.9.5.tar.gz', repos = NULL, type = "source")
# library(ggrepel)


# Index
index = read.xlsx('index.all.xlsx')
index = index |>
  select(outcome, ICD, name, Ncases) |>
  rename(BBJ = outcome)

# 添加 FEV1、FVC、FEV1_FVC的index信息
index = rbind(index, data.frame(
  BBJ = c('FEV1', 'FVC', 'FEV1_FVC'),
  ICD = rep('LungFunction', 3),
  name = c(
    'Forced expiratory volume',
    'Forced vital capacity',
    'FEV1/FVC ratio'
  ),
  Ncases = rep(20000, 3)
))

# data
mr = fread('mr_res_for_plot.csv')
snp = fread('snp_res_for_plot.csv')
# snp info index
snp_info = read.xlsx('lead_snps_160_phenotypes-v2.xlsx')
snp_info = snp_info[snp_info$SNP %in% snp$exp_name, ]
table(snp$exp_name %in% snp_info$SNP)# ALL T
snp_info = rename(snp_info, exp_name = SNP)
snp_info = snp_info |>
  filter(!(CHR == 6 & BP > 28477797 & CHR == 6 & BP < 33448354))
snp_info_not_in_HLA = unique(snp_info$exp_name)# 373

# 判断HLA
snp = snp[snp$exp_name %in% snp_info_not_in_HLA, ]
# 整理并合并数据
table(duplicated(index$BBJ))# ALL F
index_mr = index[index$BBJ %in% mr$BBJ, ]
index_snp = index[index$BBJ %in% snp$BBJ, ]

# 合并
mr = merge(mr, index_mr, by = 'BBJ', all.x = T)
snp = merge(snp, index_snp, by = 'BBJ', all.x = T)
mr$plot_group = 'mr'
snp$plot_group = 'snp'
dat = rbind(mr, snp)
unique(dat$BBJ)
unique(mr$BBJ)
unique(snp$BBJ)

# ------- 绘图 -------
dat = dat[dat$Ncases > 500,]
unique(dat$ICD)
dat$ICD_group = gsub("[^a-zA-Z]", "", map_chr(dat$ICD, ~ paste(unique(str_split(., "")[[1]]), collapse = "")))
dat$ICD_group <- ifelse(dat$ICD_group %in% c("A", "B"), "A/B", dat$ICD_group)
dat$ICD_group <- ifelse(dat$ICD == "D25", "C", dat$ICD_group)
unique(dat$ICD_group)
dat$ICD_group <- factor(dat$ICD_group,
                        levels = c(
                          "A/B", "C", "D",
                          "E", "F", "G",
                          "H", "I", "J",
                          "K", "L", "M",
                          "N", "O",
                          "Q", "R", "S",
                          "Z", 'LungFctio'
                        ),
                        labels = c(
                          "A00-B99: Infectious/parasitic",
                          "C00-D48: Neoplasms",
                          "D50-D89: Blood/blood-forming",
                          "E: Endocrine/nutritional/metabolic",
                          "F: Mental/behavioural",
                          "G: Nervous system",
                          "H: Eye/ear and adnexa",
                          "I: Circulatory system",
                          "J: Respiratory system",
                          "K: Digestive system",
                          "L: Skin/subcutaneous tissue",
                          "M: Musculoskeletal/tissue",
                          "N: Genitourinary system",
                          "O: Pregnancy/childbirth/puerperium",
                          "Q: Malformations/chromosoma",
                          "R: Other",
                          "S: Injury/poisoning",
                          "Z: Health services",
                          'Lung Function'
                        )
)
unique(dat[dat$plot_group == 'snp','exp_name'])# 373

# ----- draw plot -----
library(ggrepel)
library(ggbreak)
library(scales)

data = dat
data = data[!is.na(data$P),]
data$P = -log10(data$P)
data[data$plot_group == 'mr', 'P'] = 0-data[data$plot_group == 'mr', 'P']
# 绘图
data <- data[order(data$ICD_group, decreasing = T), ]
data <- data %>%
  group_by(ICD_group) %>%
  arrange(ICD_group, BBJ)
out_min <- data %>%
  group_by(BBJ) %>%
  reframe(
    max_p = max(P),
    min_p = min(P)
  )
out <- data[, c("BBJ", "ICD_group", 'name')]
out <- unique(out)
out_min <- merge(out, out_min, by = "BBJ", all.x = T)
out_min <- out_min[order(out_min$ICD_group, decreasing = F), ]
out_min <- out_min %>%
  arrange(ICD_group, min_p)

# out_min sig level cut
#  -log10(0.05/308)=3.789581，SNP
#  -log10(0.05/99)=3.296665，MR
data = data[data$P < 30,]
label_note = out_min[out_min$max_p > 3.789581 | out_min$min_p < -3.296665 ,]
data$BBJ <- factor(data$BBJ, levels = out_min$BBJ)
label_note_mr = label_note
label_note_snp = label_note
# 上方标注
label_note_snp$min_p = NULL
label_note_snp =  label_note_snp[label_note_snp$max_p > 3.789581,]
label_note_snp = label_note_snp[which(label_note_snp$name %in% c(
  "Colorectal cancer", "Esophageal cancer",
  "Type 2 diabetes", "Ischemic stroke", "Myocardial infarction",
   "Cholelithiasis", "Lung cancer", "Chronic obstructive pulmonary disease",
  "Interstitial lung disease", "Asthma", "Forced expiratory volume",
  "Forced vital capacity", "FEV1/FVC ratio"
)),]
# 下方标注
label_note_mr$max_p = NULL
label_note_mr = label_note_mr[label_note_mr$min_p <  -3.296665,]
label_note_mr = label_note_mr[which(label_note_mr$name %in% c(
  "Type 2 diabetes", 'Depression', 'Ischemic stroke', 'Myocardial infarction',
  "Lung cancer", "Chronic obstructive pulmonary disease",
  "Interstitial lung disease", "Asthma", "Forced expiratory volume",
  "Forced vital capacity", "FEV1/FVC ratio"
)),]
snp |>
  group_by(BBJ) |>
  reframe(n = n())# 308

fwrite(snp, file = 'lead_snp_no_hla_308BBJ.csv')


# PLOT
p = ggplot(data, aes(x = BBJ, y = P, color = ICD_group)) +
  geom_point(size = 1.2, shape = 19) +
  theme_classic() +
  geom_abline(slope = 0, intercept = 3.789581, col = "#cf142b", linetype = 2, linewidth = 0.6) +
  geom_abline(slope = 0, intercept = 0, col = "black", linetype = 1) +
  geom_abline(slope = 0, intercept = -3.296665, col = "#cf142b", linetype = 2, linewidth = 0.6) +
  scale_color_manual(values = c(
    "#eebca4",  "#f08d7b",  "#f06559",  "#e64537",  "#ca2829",  "#aa1524", 
             "#619cc5",  "#3076ad",  "#0a508f", 
             "#FFD92F","#ffbf00","#ecbe34","#F7B778", "#eca62a",
             "#c7dfc3", "#a2d099", "#59a861","#4DAF4A", "#237339" ,"#008080","#145a30"
  )) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.9, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14), # 调整图例项间距
    legend.margin = margin(t = 0), # 图例与绘图区域的上方间距
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line.y.right = element_blank(), # 禁用右边的 Y 轴线
    axis.ticks.y.right = element_blank(), # 禁用右边的 Y 轴刻度
    axis.text.y.right = element_blank()   # 禁用右边的 Y 轴标签
  ) +
  # scale_x_discrete(expand=c(0.5, 0.01))+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "Phenotype", y = "-log10(P)", color = "Phenotypic category") +
  scale_y_break(c(-40, -16)) +  # 添加第一个断裂
  scale_y_break(c(-15, -6)) +     # 添加第二个断裂
  scale_y_break(c(8.5, 12.5)) +     # 添加第三个断裂
  scale_y_continuous(breaks = c(-40, -15, -5, 0, 5, 15, 20), limits = c(-42, 20))+
  geom_text_repel(
    data = unique(subset(data, plot_group == 'mr' & P %in% label_note_mr$min_p)),
    nudge_x = 0, nudge_y = -0.2, angle = 0, min.segment.length = 0.1,
    colour = "#000000",
    force = 1,
    size = 4.5, # 调整文本大小
    aes(label = name)
  )+
  geom_text_repel(
    data = unique(subset(data, plot_group == 'snp' & P %in% label_note_snp$max_p)),
    nudge_x = 1, nudge_y = 0.3, angle = 0, min.segment.length = 0.5,
    colour = "#000000",
    size = 4.5, # 调整文本大小
    aes(label = name)
  )
ggsave(p,file = "BBJ_MR_SNP_plot_new.pdf",width = 20, height = 10)

###在ai里进一步调整修改###































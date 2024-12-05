rm(list = ls())

library(data.table)
library(tidyverse)

load("data.Rdata")

# 载入主成分数据
pc <- fread("UKBiobank_PCA.table") #ukb的PC，可以直接从ukb的tab数据中提取
pc <- pc[pc$f.eid %in% data$f.eid, ]
data <- merge(data, pc, by = "f.eid")

# 建模
fit1 <- as.formula(
    paste0(
        "AGER ~ age + sex + BMI + smoking_status + packyear + race + ",
        paste0("PC", c(1:10), collapse = " + ")
    )
)


model1 <- glm(fit1, data = data)
data$predict_AGER <- predict(model1)
data$AGER_residual <- data$AGER - data$predict_AGER

# model P value
fit_p <- as.formula(
  paste0(
    "AGER ~ rs41268920_A + age + sex + BMI + smoking_status + packyear + race + ",
    paste0("PC", c(1:10), collapse = " + ")
  )
)
model <- glm(fit_p, data = data)



library(ggbeeswarm)
library(ggsci)

# 绘图
# hard-call
get_hard_call <- function(data) {
    data <- data %>%
        mutate(
            rs41268920_Group = case_when(
                0.5 * abs(rs41268920_A - round(rs41268920_A)) > 0.1 ~ "Missing",
                rs41268920_A > 1.9 ~ "AA",
                rs41268920_A >= 0.9 & rs41268920_A < 1.1 ~ "AC",
                rs41268920_A < 0.1 ~ "CC",
                .default = "Missing"
            )
        )
    data$rs41268920_Group <- factor(data$rs41268920_Group, levels = c("AA", "AC", "CC", "Missing"))
    table(data$rs41268920_Group)
    data <- data[data$rs41268920_Group != "Missing", ]
    return(data)
}
data <- get_hard_call(data = data)
data_ns <- get_hard_call(data = data_ns)
data_es <- get_hard_call(data = data_es)

get_plot <- function(data) {
    p_temp <- ggplot(data, aes(AGER_residual, x = rs41268920_Group, fill = rs41268920_Group)) +
        geom_quasirandom(data = data[data$rs41268920_Group != "Missing"], aes(color = rs41268920_Group), size = 2, alpha = 0.25, varwidth = TRUE) +
        geom_boxplot(aes(color = rs41268920_Group),
            alpha = 0.25, width = 0.15,
            color = "black",
            position = position_dodge(width = 0.5)
        ) +
        stat_boxplot(geom = "errorbar", width = 0.25, position = position_dodge(width = 0.5)) +
        geom_violin(
            aes(
                x = rs41268920_Group,
                fill = rs41268920_Group
            ),
            color = "lightgray",
            alpha = 0.5, width = 0.5
        ) +
        scale_fill_npg() +
        scale_color_npg() +
        labs(
            x = "rs41268920", y = "AGER Protein expression",
            color = "rs41268920",
            fill = "rs41268920"
        ) +
        theme_classic() +
        stat_summary(fun.y = "median", geom = "line", size = 0.8, aes(group = 1), color = "gray") +
        stat_summary(fun.y = "median", geom = "point", size = 1)
    return(p_temp)
}
p <- get_plot(data)
p_ns <- get_plot(data_ns)
p_es <- get_plot(data_es)


ggsave(p,
    file = "p_rs41268920_AGER_pQTL_all.pdf",
    width = 3.5, height = 3.5
)
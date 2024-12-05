rm(list = ls())

library(tidyverse)
library(survival)
library(data.table)

load("AGER_LC_COPD_ILD_Asthma_Cox_data_all_clean.Rdata")

# 基本参数构建
data$LC_time <- as.numeric(data$LC_date - data$date_assessment)
data$COPD_time <- as.numeric(data$COPD_date - data$date_assessment)
data$ILD_time <- as.numeric(data$ILD_date - data$date_assessment)
data$Asthma_time <- as.numeric(data$Asthma_date - data$date_assessment)

data[is.na(data$COPD_base), 'COPD_base'] = 0

# 构建函数完成连续主效应分析
get_hr <- function(data, event_name, time_name) {
    # 全效应
    m_main <- coxph(
        Surv(data[[time_name]], data[[event_name]]) ~ result +
            age + sex + BMI + income + education + smoking_status +
            packyear + alcohol_status + family_cancer + race,
        data = data
    )
    res_main <- data.frame(
        Terms = paste0(event_name, ": Total"),
        Case_total = paste0(m_main$nevent, "/", m_main$n),
        Incidence =  round((m_main$nevent / (sum(data[[time_name]])/365.25))*1e05,2),
        Beta = round(summary(m_main)$coef[1, 1], 2),
        SE = round(summary(m_main)$coef[1, 3], 2),
        HR = round(summary(m_main)$conf.int[1, 1], 2),
        CI = paste0(
            round(summary(m_main)$conf.int[1, 3], 2),
            "-",
            round(summary(m_main)$conf.int[1, 4], 2)
        ),
        P = signif(summary(m_main)$coef[1, 5], digits = 3)
    )

    # 亚组数据
    data_ns <- data[data$smoking_status == "Never", ]
    data_es <- data[data$smoking_status != "Never", ]

    # 非吸烟
    m_ns <- coxph(
        Surv(data_ns[[time_name]], data_ns[[event_name]]) ~ result +
            age + sex + BMI + income + education +
            alcohol_status + family_cancer + race,
        data = data_ns
    )
    res_ns <- data.frame(
        Terms = paste0(event_name, ": Non smoker"),
        Case_total = paste0(m_ns$nevent, "/", m_ns$n),
        Incidence =  round((m_ns$nevent / (sum(data_ns[[time_name]])/365.25))*1e05,2),
        Beta = round(summary(m_ns)$coef[1, 1], 2),
        SE = round(summary(m_ns)$coef[1, 3], 2),
        HR = round(summary(m_ns)$conf.int[1, 1], 2),
        CI = paste0(
            round(summary(m_ns)$conf.int[1, 3], 2),
            "-",
            round(summary(m_ns)$conf.int[1, 4], 2)
        ),
        P = signif(summary(m_ns)$coef[1, 5], digits = 3)
    )
    # 吸烟者
    m_es <- coxph(
        Surv(data_es[[time_name]], data_es[[event_name]]) ~ result +
            age + sex + BMI + income + education + packyear +
            alcohol_status + family_cancer + race,
        data = data_es
    )
    res_es <- data.frame(
        Terms = paste0(event_name, ": Ever smoker"),
        Case_total = paste0(m_es$nevent, "/", m_es$n),
        Incidence =  round((m_es$nevent / (sum(data_es[[time_name]])/365.25))*1e05,2),
        Beta = round(summary(m_es)$coef[1, 1], 2),
        SE = round(summary(m_es)$coef[1, 3], 2),
        HR = round(summary(m_es)$conf.int[1, 1], 2),
        CI = paste0(
            round(summary(m_es)$conf.int[1, 3], 2),
            "-",
            round(summary(m_es)$conf.int[1, 4], 2)
        ),
        P = signif(summary(m_es)$coef[1, 5], digits = 3)
    )

    # 合并结果
    rbind(
        res_main, res_ns, res_es
    )
}

# 分析与合并
res <- rbind(
    get_hr(data[data$COPD_base == '0',], event_name = "LC", time_name = "LC_time"),
    get_hr(data[data$COPD_base == '0',], event_name = "COPD", time_name = "COPD_time"),
    get_hr(data[data$COPD_base == '0',], event_name = "ILD", time_name = "ILD_time"),
    get_hr(data[data$COPD_base == '0',], event_name = "Asthma", time_name = "Asthma_time")
)
res

data = data%>%
  select(
    -c(vitb6_mean, fev1, fvc, Z_fev1_fvc, vitb6_1st,
       vitb6_2nd, vitb6_3rd, vitb6_4th, vitb6_5th, cancer, fev1_fvc)
  )
# 抽样核对原始数据
data[666,]# 与原始数据核对一致 ALL T
save(data, file = 'AGER_LC_COPD_ILD_Asthma_Cox_data_all_clean.Rdata')

# 写出
# 输出到word
library(officer)
library(flextable)
save_dataframes_to_word <- function(path, ...) {
    dataframes <- list(...)
    df_titles <- lapply(dataframes, function(df) attr(df, "title"))
    names_list <- as.character(substitute(list(...))[-1])
    doc <- read_docx()
    for (i in seq_along(dataframes)) {
        title <- ifelse(is.null(df_titles[[i]]), names_list[i], df_titles[[i]])
        doc <- doc %>%
            body_add_par(value = title, style = "heading 1") %>%
            body_add_flextable(flextable(dataframes[[i]]))

        if (i < length(dataframes)) {
            doc <- doc %>%
                body_add_par(value = "") %>%
                body_add_par(value = "") %>%
                body_add_par(value = "")
        }
    }
    print(doc, target = path)
}
# 定义表格标题
attr(res, "title") <- "AGER-LC_COPD_ILD_Asthma_cox_regression_result"
# 输出
save_dataframes_to_word(
    path = "AGER-LC_COPD_ILD_Asthma_cox_regression_result.docx",
    res
)

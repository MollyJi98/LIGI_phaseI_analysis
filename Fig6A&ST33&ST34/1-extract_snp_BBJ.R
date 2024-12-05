R
rm(list = ls())

library(tidyverse)
library(data.table)
library(readxl)

snp <- readxl::read_xlsx("lead_snps_160_phenotypes-v2.xlsx")
snp <- unique(snp)

files <- fread("valid_folders_name.txt", header = F) #BBJ表型
colnames(files) <- "folder_name"
files$path <- paste0("./BBJ_ALL_END_SUM_STAT/uncompressed/", files$folder_name)
for (i in 1:nrow(files)) {
  dir_command <- paste0("ls ", files$path[i], "| grep .auto.txt")
  files$gzfile_name[i] <- system(dir_command, intern = T)
}
files$gzpath <- paste0(files$path, "/", files$gzfile_name)
files$name <- str_extract(files$folder_name, "(?<=BBJ\\.).*(?=\\.v1)")

# 声明
dat_com <- data.frame(
  rsID = snp$SNP,
  CHR_BP = snp$CHR_BP
)
# loop
for (i in 1:nrow(files)) {
  dat <- fread(files$gzpath[i])

  dat$CHR_BP <- paste0(dat$CHR, ":", dat$POS)
  dat <- dat[dat$CHR_BP %in% dat_com$CHR_BP, ]

  if (nrow(dat) == 0) {
    dat <- data.frame(
      CHR_BP = dat_com$CHR_BP
    )
    cat("---", files$name[i], " non snp been found ---- \n")
  }

  if (nrow(dat) > 0) {
    A1_col <- paste0(files$name[i], "_A1")
    A2_col <- paste0(files$name[i], "_A2")
    A1_freq_col <- paste0(files$name[i], "_A1_freq")
    Beta_col <- paste0(files$name[i], "_Beta")
    P_col <- paste0(files$name[i], "_Pval")
    dat <- dat %>%
      select(
        CHR_BP, Allele2, Allele1, AF_Allele2, BETA, p.value
      ) %>%
      rename(
        CHR_BP = CHR_BP,
        !!sym(A1_col) := Allele2,
        !!sym(A2_col) := Allele1,
        !!sym(A1_freq_col) := AF_Allele2,
        !!sym(Beta_col) := BETA,
        !!sym(P_col) := p.value,
      )
    cat("---", files$name[i], " has been processed! ---- \n")
  }
  dat_com <- merge(dat_com, dat, by = "CHR_BP", all.x = T)
}
write.table(dat_com, file = "lean_snp_matrix.txt", quote = F, row.names = F)

# 本地
rm(list = ls())

library(tidyverse)
library(data.table)

dat <- fread("lean_snp_matrix.txt")
table(is.na(dat$PT_A1)) # 318

dat_A1 <- dat %>%
  select(rsID, ends_with("_A1")) %>%
  pivot_longer(
    cols = ends_with("_A1"),
    names_to = "variable",
    values_to = "A1"
  ) %>%
  mutate(
    variable = str_replace(variable, "_A1", ""),
    merID = paste0(rsID, "_", variable)
  ) %>%
  rename(
    outcome = variable
  )
dat_A2 <- dat %>%
  select(rsID, ends_with("_A2")) %>%
  pivot_longer(
    cols = ends_with("_A2"),
    names_to = "variable",
    values_to = "A2"
  ) %>%
  mutate(
    variable = str_replace(variable, "_A2", ""),
    merID = paste0(rsID, "_", variable)
  ) %>%
  select(
    merID, A2
  )
dat_A1_freq <- dat %>%
  select(rsID, ends_with("_freq")) %>%
  pivot_longer(
    cols = ends_with("_freq"),
    names_to = "variable",
    values_to = "A1_freq"
  ) %>%
  mutate(
    variable = str_replace(variable, "_A1_freq", ""),
    merID = paste0(rsID, "_", variable)
  ) %>%
  select(
    merID, A1_freq
  )
dat_Beta <- dat %>%
  select(rsID, ends_with("_Beta")) %>%
  pivot_longer(
    cols = ends_with("_Beta"),
    names_to = "variable",
    values_to = "Beta"
  ) %>%
  mutate(
    variable = str_replace(variable, "_Beta", ""),
    merID = paste0(rsID, "_", variable)
  ) %>%
  select(
    merID, Beta
  )
dat_Pval <- dat %>%
  select(rsID, ends_with("_Pval")) %>%
  pivot_longer(
    cols = ends_with("_Pval"),
    names_to = "variable",
    values_to = "Pval"
  ) %>%
  mutate(
    variable = str_replace(variable, "_Pval", ""),
    merID = paste0(rsID, "_", variable)
  ) %>%
  select(
    merID, Pval
  )

data <- merge(dat_A1, dat_A2, by = "merID", all.x = T)
data <- merge(data, dat_A1_freq, by = "merID", all.x = T)
data <- merge(data, dat_Beta, by = "merID", all.x = T)
data <- merge(data, dat_Pval, by = "merID", all.x = T)
data <- data[!is.na(data$Beta), ]
unique(data$rsID) # 318
data$merID <- NULL

library(readxl)
index <- read_xlsx("index.all.xlsx")
index$outcome <- str_extract(index$BBJ, "(?<=BBJ\\.).*(?=\\.v1)")
index <- index[, c("Ncases", "Ncontrols", "outcome", "CN")]
data <- merge(data, index, by = "outcome", all.x = T)
library(openxlsx)
write.xlsx(data, file = "lead_snp_BBJ_info_all_renew.xlsx")

############在excel里整理结果#################################################################
############最终报导的位点是提取到的318个位点-MHC区域的10个位点，即308个lead snp的效应情况############
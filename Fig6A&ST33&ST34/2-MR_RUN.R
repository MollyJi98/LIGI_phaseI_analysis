# 准备结局文件索引
rm(list = ls())
library(data.table)
library(tidyverse)
# 以结局为主，构建结局文件固定索引，然后分10进程并行完成

end_files <- data.frame(
    end_folder = system("ls ./data/for_MR", intern = T)
) %>%
    mutate(
        end_name = str_extract(end_folder, "(?<=BBJ\\.).*(?=\\.v1)"),
        end_folder_path = paste0("./data/for_MR/", end_folder)
    )
seq(1, 159, by = 16) # 1  17  33  49  65  81  97 113 129 145
end_files_01 <- end_files[1:17, ]
end_files_02 <- end_files[18:33, ]
end_files_03 <- end_files[34:49, ]
end_files_04 <- end_files[50:65, ]
end_files_05 <- end_files[66:81, ]
end_files_06 <- end_files[82:97, ]
end_files_07 <- end_files[98:113, ]
end_files_08 <- end_files[114:129, ]
end_files_09 <- end_files[130:145, ]
end_files_10 <- end_files[146:159, ]
end_files_all <- end_files
# 写出
write.table(end_files_all, file = "end_files_all.txt", quote = F, row.names = F)
write.table(end_files_01, file = "end_files_01.txt", quote = F, row.names = F)
write.table(end_files_02, file = "end_files_02.txt", quote = F, row.names = F)
write.table(end_files_03, file = "end_files_03.txt", quote = F, row.names = F)
write.table(end_files_04, file = "end_files_04.txt", quote = F, row.names = F)
write.table(end_files_05, file = "end_files_05.txt", quote = F, row.names = F)
write.table(end_files_06, file = "end_files_06.txt", quote = F, row.names = F)
write.table(end_files_07, file = "end_files_07.txt", quote = F, row.names = F)
write.table(end_files_08, file = "end_files_08.txt", quote = F, row.names = F)
write.table(end_files_09, file = "end_files_09.txt", quote = F, row.names = F)
write.table(end_files_10, file = "end_files_10.txt", quote = F, row.names = F)

# RUN MR
# screen -r 208424.JC_BBJ_MR_RUN_10
# screen -r 205372.JC_BBJ_MR_RUN_09
# screen -r 202142.JC_BBJ_MR_RUN_08
# screen -r 198854.JC_BBJ_MR_RUN_07
# screen -r 196036.JC_BBJ_MR_RUN_06
# screen -r 192866.JC_BBJ_MR_RUN_05
# screen -r 188357.JC_BBJ_MR_RUN_04
# screen -r 183872.JC_BBJ_MR_RUN_03
# screen -r 177627.JC_BBJ_MR_RUN_02
# screen -r 170199.JC_BBJ_MR_RUN_01
rm(list = ls())
library(data.table)
library(tidyverse)
library(TwoSampleMR)

end_files <- fread("end_files_01.txt")
# end_files <- fread("end_files_02.txt")
# end_files <- fread("end_files_03.txt")
# end_files <- fread("end_files_04.txt")
# end_files <- fread("end_files_05.txt")
# end_files <- fread("end_files_06.txt")
# end_files <- fread("end_files_07.txt")
# end_files <- fread("end_files_08.txt")
# end_files <- fread("end_files_09.txt")
# end_files <- fread("end_files_10.txt")

# RUN
for (i in 1:nrow(end_files)) {
    # 构建暴露和结局数据索引
    mr_path <- end_files$end_folder_path[i]
    exp_name_com <- paste0("ls ", end_files$end_folder_path[i], "| grep .for.mr.exp")
    exp_name_com <- system(exp_name_com, intern = T)
    exp_name_com <- str_extract(exp_name_com, ".*(?=\\.for\\.mr\\.exp)")
    mr_files <- data.frame(
        mr_path = mr_path,
        end_name = exp_name_com
    ) %>%
        mutate(
            mr_exp_path = paste0(mr_path, "/", end_name, ".for.mr.exp"),
            mr_out_path = paste0(mr_path, "/", end_name, ".for.mr.out")
        )

    # 声明
    res_com <- NULL
    het_com <- NULL
    ple_com <- NULL
    # MR loop
    for (j in 1:nrow(mr_files)) {
        # read data
        exp_dat <- read_exposure_data(
            mr_files$mr_exp_path[j],
            snp_col = "rsID",
            effect_allele_col = "A1",
            other_allele_col = "A2",
            eaf_col = "A1_freq",
            beta_col = "Beta",
            se_col = "SE",
            pval_col = "P"
        )
        out_dat <- read_outcome_data(
            mr_files$mr_out_path[j],
            snp_col = "rsID",
            effect_allele_col = "A1",
            other_allele_col = "A2",
            eaf_col = "A1_freq",
            beta_col = "Beta",
            se_col = "SE",
            pval_col = "P"
        )
        # pheno
        exp_dat$exposure <- mr_files$end_name[j]
        out_dat$outcome <- end_files$end_name[i]
        # harmonise
        dat <- harmonise_data(exp_dat, out_dat)

        # 5e-08
        dat <- dat[dat$pval.exposure < 5e-08, ]

        # mr
        res <- mr(dat)
        # het
        het <- mr_heterogeneity(dat)
        # ple
        ple <- mr_pleiotropy_test(dat)

        # 结果合并
        res_com <- rbind(res_com, res)
        het_com <- rbind(het_com, het)
        ple_com <- rbind(ple_com, ple)

        # 写出mr run data
        dat_output_dir <- paste0(end_files$end_folder[i], "/", mr_files$end_name[j], ".mr.run.dat")
        write.table(dat, file = dat_output_dir, quote = F, row.names = F)

        # cat
        cat("---- ", mr_files$end_name[j], " has been processed! ---- \n")
    }
    dir_res <- paste0(end_files$end_name[i], ".mr.res.csv")
    dir_het <- paste0(end_files$end_name[i], ".mr.het.csv")
    dir_ple <- paste0(end_files$end_name[i], ".mr.ple.csv")
    write.csv(res_com, file = dir_res)
    write.csv(res_com, file = dir_het)
    write.csv(res_com, file = dir_ple)
    # cat
    cat("======== ", end_files$end_name[i], " has been done! ====== \n")
}

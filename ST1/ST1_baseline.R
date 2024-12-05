#################baseline#################
R
library(data.table)
library(tableone)

#load file
cov=fread('COV_ZJ_JS_imaging_39334.txt')
cov=cov[,-c(1:2)]
sample=fread('ZJ2.8_JS_genotype_merge_for_BOLTLMM.fam')
head(sample)
sample <- sample[-1,]

table(sample$V1 %in% cov$FID)
cov1 <- subset(cov,cov$FID %in% sample$V1)
table(is.na(cov1$sex)) #35469

cov2 <- subset(cov1,is.na(cov1$sex)==F)
str(cov2)
cov2$region <- as.factor(cov2$region)
cov2$sex <- as.factor(cov2$sex)
cov2$Education <- as.factor(cov2$Education)
cov2$Occupational_hazard <- as.factor(cov2$Occupational_hazard)
cov2$smoke <- as.factor(cov2$smoke)
cov2$smoke2 <- as.factor(ifelse(cov2$smoke==0,0,ifelse(cov2$smoke==1|cov2$smoke==2,1,999)))
cov2$ETS_at_home <- as.factor(cov2$ETS_at_home)
cov2$ETS_at_work <- as.factor(cov2$ETS_at_work)
cov2$Exercise <- as.factor(cov2$Exercise)
cov2$Alcohol_drinking <- as.factor(cov2$Alcohol_drinking)
cov2$family_history_LC <- as.factor(cov2$family_history_LC)
cov2$family_history_LC_firstdegree <- as.factor(cov2$family_history_LC_firstdegree)
cov2$cancer_history <- as.factor(cov2$cancer_history)

var <- names(cov2)[c(4:10,12:24)]
result1 <- CreateTableOne(cov2, vars = var, strata = 'sex')
tab4Mat <- print(result1,printToggle = FALSE,catDigits = 2)

result2 <- CreateTableOne(data = cov2, vars = var)
tab4Mat2 <- print(result2,printToggle = FALSE,catDigits = 2)

tab4Mat3 <- cbind(tab4Mat,tab4Mat2)
write.table(tab4Mat3,'baseline.xlsx',quote = F,sep = '\t',row.names = T)



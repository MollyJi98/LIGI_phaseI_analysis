rm(list=ls())

library(readxl)
library(forestploter)
library(ggplot2)
library(grid)
library(ggsci)
library(openxlsx) #4.2.3
library(forestplot) #1.10.1


###########################rs41268920###########################
rs41268920 <- read_excel("plot_file.xlsx",sheet="rs41268920")
desired_order <- c("NSCLC","COPD","ILD","ASTHMA","FEV1","FVC","FEV1/FVC")
rs41268920$outcome <- factor(rs41268920$outcome, levels = rev(desired_order))

##4种疾病
rs41268920_diseasse <- rs41268920[1:4,]
tabletext <- cbind(c("outcome",as.character(rs41268920_diseasse$outcome)),
                   c("Hazard Ratio (95% CI)",as.character(rs41268920_diseasse$`Hazard Ratio (95% CI)`)),
                   c("P value",as.character(rs41268920_diseasse$`P value`)))

OR2<-c(NA,rs41268920_diseasse$hr)
IL2<-c(NA,rs41268920_diseasse$low)
UL2<-c(NA,rs41268920_diseasse$up)
dat3=as.data.frame(cbind(OR2,IL2,UL2))


pdf(file = "rs41268920_diseasse.pdf",width = 8, height = 8,bg="transparent")
forestplot(labeltext = tabletext[,c(1:3)],
           dat3$OR2,dat3$IL2,dat3$UL2,
           graph.pos=3,
           boxsize   = 0.4,
           line.margin=1.0,
           lineheight=unit(6.50,'mm'),
           colgap=unit(1.00,'cm'),
           zero      = 1,cex=2,
           xlog      = F,
           clip=c(0.8,1.2),
           txt_gp=fpTxtGp(ticks=gpar(cex=0.5),title =gpar(cex=0.5),xlab=gpar(cex=0.01)),
           xticks=(c(0.8,0.9,1.00,1.1,1.2)),
           col = fpColors(lines=c("#747474", "#000000","#747474", "#000000"), box=c("#747474", "#000000","#747474", "#000000"),zero = "#000000"),
           graphwidth=unit(3.5,'cm'),
           zero.line = gpar(lty = "dashed"),
           refline_lwd = 1, 
           refline_lty ="dashed",vertline_lwd = 1, vertline_lty ="dashed",
)

dev.off()

##3种肺功能
rs41268920_diseasse <- rs41268920[5:7,]
tabletext <- cbind(c("outcome",as.character(rs41268920_diseasse$outcome)),
                   c("Hazard Ratio (95% CI)",as.character(rs41268920_diseasse$`Hazard Ratio (95% CI)`)),
                   c("P value",as.character(rs41268920_diseasse$`P value`)))

OR2<-c(NA,rs41268920_diseasse$hr)
IL2<-c(NA,rs41268920_diseasse$low)
UL2<-c(NA,rs41268920_diseasse$up)
dat3=as.data.frame(cbind(OR2,IL2,UL2))


pdf(file = "rs41268920_lungfunction.pdf",width = 8, height = 8,bg="transparent")
forestplot(labeltext = tabletext[,c(1:3)],
           dat3$OR2,dat3$IL2,dat3$UL2,
           graph.pos=3,
           boxsize   = 0.4,
           line.margin=1.0,
           lineheight=unit(6.50,'mm'),
           colgap=unit(1.00,'cm'),
           zero      = 0,cex=2,
           xlog      = F,
           clip=c(-0.06,0.05),
           txt_gp=fpTxtGp(ticks=gpar(cex=0.5),title =gpar(cex=0.5),xlab=gpar(cex=0.01)),
           xticks=(c(-0.05,-0.025,0,0.025,0.05)),
           col = fpColors(lines=c("#747474", "#000000","#747474", "#000000"), box=c("#747474", "#000000","#747474", "#000000"),zero = "#000000"),
           graphwidth=unit(3.5,'cm'),
           zero.line = gpar(lty = "dashed"),
           refline_lwd = 1, 
           refline_lty ="dashed",vertline_lwd = 1, vertline_lty ="dashed",
)

dev.off()







###########################rs4693974###########################
rs4693974 <- read_excel("准备文件2.xlsx",sheet="rs4693974")
desired_order <- c("NSCLC","COPD","ILD","ASTHMA","FEV1","FVC","FEV1/FVC")
rs4693974$outcome <- factor(rs4693974$outcome, levels = rev(desired_order))

##4种疾病
rs4693974_diseasse <- rs4693974[1:4,]
tabletext <- cbind(c("outcome",as.character(rs4693974_diseasse$outcome)),
                   c("Hazard Ratio (95% CI)",as.character(rs4693974_diseasse$`Hazard Ratio (95% CI)`)),
                   c("P value",as.character(rs4693974_diseasse$`P value`)))

OR2<-c(NA,rs4693974_diseasse$hr)
IL2<-c(NA,rs4693974_diseasse$low)
UL2<-c(NA,rs4693974_diseasse$up)
dat3=as.data.frame(cbind(OR2,IL2,UL2))


pdf(file = "rs4693974_diseasse.pdf",width = 8, height = 8,bg="transparent")
forestplot(labeltext = tabletext[,c(1:3)],
           dat3$OR2,dat3$IL2,dat3$UL2,
           graph.pos=3,
           boxsize   = 0.4,
           line.margin=1.0,
           lineheight=unit(6.50,'mm'),
           colgap=unit(1.00,'cm'),
           zero      = 1,cex=2,
           xlog      = F,
           clip=c(0.75,1.2),
           txt_gp=fpTxtGp(ticks=gpar(cex=0.5),title =gpar(cex=0.5),xlab=gpar(cex=0.01)),
           xticks=(c(0.8,0.9,1.00,1.1,1.2)),
           col = fpColors(lines=c("#747474", "#000000","#747474", "#000000"), box=c("#747474", "#000000","#747474", "#000000"),zero = "#000000"),
           graphwidth=unit(3.5,'cm'),
           zero.line = gpar(lty = "dashed"),
           refline_lwd = 1, 
           refline_lty ="dashed",vertline_lwd = 1, vertline_lty ="dashed",
)

dev.off()

##3种肺功能
rs4693974_diseasse <- rs4693974[5:7,]
tabletext <- cbind(c("outcome",as.character(rs4693974_diseasse$outcome)),
                   c("Hazard Ratio (95% CI)",as.character(rs4693974_diseasse$`Hazard Ratio (95% CI)`)),
                   c("P value",as.character(rs4693974_diseasse$`P value`)))

OR2<-c(NA,rs4693974_diseasse$hr)
IL2<-c(NA,rs4693974_diseasse$low)
UL2<-c(NA,rs4693974_diseasse$up)
dat3=as.data.frame(cbind(OR2,IL2,UL2))


pdf(file = "rs4693974_lungfunction.pdf",width = 8, height = 8,bg="transparent")
forestplot(labeltext = tabletext[,c(1:3)],
           dat3$OR2,dat3$IL2,dat3$UL2,
           graph.pos=3,
           boxsize   = 0.4,
           line.margin=1.0,
           lineheight=unit(6.50,'mm'),
           colgap=unit(1.00,'cm'),
           zero      = 0,cex=2,
           xlog      = F,
           clip=c(-0.06,0.05),
           txt_gp=fpTxtGp(ticks=gpar(cex=0.5),title =gpar(cex=0.5),xlab=gpar(cex=0.01)),
           xticks=(c(-0.05,-0.025,0,0.025,0.05)),
           col = fpColors(lines=c("#747474", "#000000","#747474", "#000000"), box=c("#747474", "#000000","#747474", "#000000"),zero = "#000000"),
           graphwidth=unit(3.5,'cm'),
           zero.line = gpar(lty = "dashed"),
           refline_lwd = 1, 
           refline_lty ="dashed",vertline_lwd = 1, vertline_lty ="dashed",
)

dev.off()

#############在ai和ppt里修改和拼图#############



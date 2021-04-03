setwd("/Users/cnhirsch/Documents/Research/Nathan/Christine_TE_methods/Analysis/LD_analysis_revisions/Dprime/")

LD2<-read.table(file="plink_results_SNPs-highest-LD-TE_Dprime.ld", head=T, sep="\t")
LD2$DP <- as.numeric(as.character(LD2$DP))
LD2$dist_to_te <- as.numeric(as.character(LD2$dist_to_te))

pdf.options(family="Helvetica")
pdf(file="TE_LDprime.pdf", pointsize=20)
bins=seq(0,1,by=0.005)
hist(LD2$DP, breaks=bins, col="gray", border="gray", axes=F, xlab="LD (D')", ylab="Frequency", main="", xlim=c(0,1), ylim=c(0,280000))
axis(1, xaxp=c(0,1,5), pos=0)
axis(2, yaxp=c(0,280000,4), pos=0, las=1)
dev.off()

pdf.options(family="Helvetica")
pdf(file="TE_distance.pdf", pointsize=20)
LD3<-subset(LD2, R2>0.9)
bins=seq(0,1000000,by=2000)
hist(LD3$dist_to_te, breaks=bins, col="gray", border="gray", axes=F, xlab="Distance (Kb)", ylab="Frequency", main="", xlim=c(0,1000000), ylim=c(0,16000))
axis(1, at=c(0, 200000, 400000, 600000, 800000, 1000000), labels=c(0, 200, 400, 600, 800, 1000), pos=0)
axis(2, yaxp=c(0,16000,4), pos=0, las=1)
dev.off()

summary(LD3)
##########################################################################################################################################################

helitron_LD<-read.table(file="TE_SNP_LDprime_helitron.txt", head=T, sep="\t")
LTR_LD<-read.table(file="TE_SNP_LDprime_LTR.txt", head=T, sep="\t")
TIR_LD<-read.table(file="TE_SNP_LDprime_TIR.txt", head=T, sep="\t")

helitron_LD$DP <- as.numeric(as.character(helitron_LD$DP))
helitron_low<-subset(helitron_LD, DP<0.5)
helitron_med<-subset(helitron_LD, DP>=0.5 & DP<=0.9)
helitron_high<-subset(helitron_LD, DP>0.9)

LTR_LD$DP <- as.numeric(as.character(LTR_LD$DP))
LTR_low<-subset(LTR_LD, DP<0.5)
LTR_med<-subset(LTR_LD, DP>=0.5 & DP<=0.9)
LTR_high<-subset(LTR_LD, DP>0.9)

TIR_LD$DP <- as.numeric(as.character(TIR_LD$DP))
TIR_low<-subset(TIR_LD, DP<0.5)
TIR_med<-subset(TIR_LD, DP>=0.5 & DP<=0.9)
TIR_high<-subset(TIR_LD, DP>0.9)

order<-c(rep("LTR", 3), rep("Helitron", 3), rep("TIR", 3))
ld<-c(rep(c("LD <0.5", "LD 0.5-0.9", "LD >0.9"), 3))
value<-c(13, 89, 114013, 4, 34, 10421, 14, 203, 118306)
data<-data.frame(order, ld, value)

data$ld<-factor(data$ld, levels=c("LD >0.9", "LD 0.5-0.9", "LD <0.5"))
library(ggplot2)


pdf.options(family="Helvetica")
pdf(file="LDprime_by_order.pdf", pointsize=30)

ggplot(data, aes(fill=ld, y=value, x=order)) + geom_bar(position="fill", stat="identity") + theme_classic() + theme(legend.title=element_blank()) + scale_fill_manual(values=c("#008dce", "#8e2043", "#077187")) + ylab("Proportion of TEs") + xlab("") + theme(text = element_text(size=30))

dev.off()


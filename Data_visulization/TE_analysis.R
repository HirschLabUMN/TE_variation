#TE-SNP LD analysis

setwd("/Users/cnhirsch/Documents/Research/Nathan/Christine_TE_methods/Analysis/LD_analysis_revisions/R2/")

LD2<-read.table(file="TE_SNP_LD_with_freq.txt", head=T, sep="\t")
LD2$R2 <- as.numeric(as.character(LD2$R2))
LD2$dist_to_te <- as.numeric(as.character(LD2$dist_to_te))

pdf.options(family="Helvetica")
pdf(file="TE_LD.pdf", pointsize=20)
bins=seq(0,1,by=0.01)
hist(LD2$R2, breaks=bins, col="gray", border="gray", axes=F, xlab="LD (r2)", ylab="Frequency", main="", xlim=c(0,1), ylim=c(0,30000))
axis(1, xaxp=c(0,1,5), pos=0)
axis(2, yaxp=c(0,30000,5), pos=0, las=1)
dev.off()

pdf.options(family="Helvetica")
pdf(file="TE_dist_all.pdf", pointsize=20)
bins=seq(0,1000000,by=2000)
hist(LD2$dist_to_te, breaks=bins, col="gray", border="gray", axes=F, xlab="Distance (Kb)", ylab="Frequency", main="", xlim=c(0,1000000), ylim=c(0,5000))
axis(1, at=c(0, 200000, 400000, 600000, 800000, 1000000), labels=c(0, 200, 400, 600, 800, 1000), pos=0)
axis(2, yaxp=c(0,5000,5), pos=0, las=1)
dev.off()

LD3<-subset(LD2, R2>0.9)
pdf.options(family="Helvetica")
pdf(file="TE_dist_r2_above_0.9.pdf", pointsize=20)
bins=seq(0,1000000,by=2000)
hist(LD3$dist_to_te, breaks=bins, col="gray", border="gray", axes=F, xlab="Distance (Kb)", ylab="Frequency", main="", xlim=c(0,1000000), ylim=c(0,4000))
axis(1, at=c(0, 200000, 400000, 600000, 800000, 1000000), labels=c(0, 200, 400, 600, 800, 1000), pos=0)
axis(2, yaxp=c(0,4000,4), pos=0, las=1)
dev.off()

summary(LD3)
nrow(subset(LD3, dist_to_te < 200000))
nrow(LD3)

nrow(subset(LD2, R2 < 0.5))
nrow(LD2)


######

helitron_LD<-read.table(file="TE_SNP_LD_helitron_with_freq.txt", head=T, sep="\t")
LTR_LD<-read.table(file="TE_SNP_LD_LTR_with_freq.txt", head=T, sep="\t")
TIR_LD<-read.table(file="TE_SNP_LD_TIR_with_freq.txt", head=T, sep="\t")

helitron_LD$R2 <- as.numeric(as.character(helitron_LD$R2))
helitron_low<-subset(helitron_LD, R2<0.5)
helitron_med<-subset(helitron_LD, R2>=0.5 & R2<=0.9)
helitron_high<-subset(helitron_LD, R2>0.9)

LTR_LD$R2 <- as.numeric(as.character(LTR_LD$R2))
LTR_low<-subset(LTR_LD, R2<0.5)
LTR_med<-subset(LTR_LD, R2>=0.5 & R2<=0.9)
LTR_high<-subset(LTR_LD, R2>0.9)

TIR_LD$R2 <- as.numeric(as.character(TIR_LD$R2))
TIR_low<-subset(TIR_LD, R2<0.5)
TIR_med<-subset(TIR_LD, R2>=0.5 & R2<=0.9)
TIR_high<-subset(TIR_LD, R2>0.9)

order<-c(rep("LTR", 3), rep("Helitron", 3), rep("TIR", 3))
ld<-c(rep(c("LD <0.5", "LD 0.5-0.9", "LD >0.9"), 3))
#value<-c(10111, 15019, 67872, 1874, 1880, 5941, 24054, 25858, 57353)
value<-c(17660, 29367, 67082, 3650, 4714, 5694, 27821, 36270, 54408)
data<-data.frame(order, ld, value)

data$ld<-factor(data$ld, levels=c("LD >0.9", "LD 0.5-0.9", "LD <0.5"))

library(ggplot2)

pdf.options(family="Helvetica")
pdf(file="LD_by_order.pdf", pointsize=30)

ggplot(data, aes(fill=ld, y=value, x=order)) + geom_bar(position="fill", stat="identity") + theme_classic() + theme(legend.title=element_blank()) + scale_fill_manual(values=c("#008dce", "#8e2043", "#077187")) + ylab("Proportion of TEs") + xlab("") + theme(text = element_text(size=30))

dev.off()


######


pdf("freq_density_helitron.pdf", pointsize=20)
plot(density(helitron_low$prop_prsent), col="#077187", xlab="Population Frequency", main="Helitron", lwd=2, yaxt="n")
lines(density(helitron_med$prop_prsent), col="#8e2043", lwd=2)
lines(density(helitron_high$prop_prsent), col="#008dce", lwd=2)
axis(2, at=c(0,2,4,6,8), labels=c(0,2,4,6,8), las=1)
legend(0, 8, c("LD >0.9", "LD 0.5-0.9", "LD <0.5"), col=c("#008dce", "#8e2043", "#077187"), lty=1, lwd=2)
dev.off()

pdf("freq_density_TIR.pdf", pointsize=20)
plot(density(TIR_low$prop_prsent), col="#077187", xlab="Population Frequency", main="TIR", lwd=2, yaxt="n")
lines(density(TIR_med$prop_prsent), col="#8e2043", lwd=2)
lines(density(TIR_high$prop_prsent), col="#008dce", lwd=2)
axis(2, at=c(0,3,6,9,12), labels=c(0,3,6,9,12), las=1)
legend(0, 13, c("LD >0.9", "LD 0.5-0.9", "LD <0.5"), col=c("#008dce", "#8e2043", "#077187"), lty=1, lwd=2)
dev.off()

pdf("freq_density_LTR.pdf", pointsize=20)
plot(density(LTR_low$prop_prsent), col="#077187", xlab="Population Frequency", main="LTR", lwd=2, yaxt="n")
lines(density(LTR_med$prop_prsent), col="#8e2043", lwd=2)
lines(density(LTR_high$prop_prsent), col="#008dce", lwd=2)
axis(2, at=c(0,4,8,12,16), labels=c(0,4,8,12,16), las=1)
legend(0, 17, c("LD >0.9", "LD 0.5-0.9", "LD <0.5"), col=c("#008dce", "#8e2043", "#077187"), lty=1, lwd=2)
dev.off()


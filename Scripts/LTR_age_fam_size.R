setwd("~/Desktop/TE paper/")
# LTR age 
WiDiv508_refB73_TEmetadata_age <- read.csv("File_S5_fixed.txt",sep = "\t") 
B73_LTRs <- subset(WiDiv508_refB73_TEmetadata_age, order == "LTR")

corr_LTR_age <- cor.test(x=B73_LTRs$fam_size, y=B73_LTRs$LTR_age, method = 'spearman')
corr_LTR_age


p <- ggplot(B73_LTRs, aes(prop_present,LTR_age))
new_color = c("#8000FFFF","#00FFFFFF","#80FF00FF","#FF0000FF")
p1 <- p + geom_bin2d(bins=150) +scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 40, by=10), colours=new_color) + 
  #stat_smooth(aes(colour="",method="losse",se = FALSE))+ theme_bw()
  stat_smooth(method="loess",se = FALSE)+ theme_bw()
              



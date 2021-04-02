setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
library(tidyr)
library(tidyverse)
library(lemon)
# heatmap 
WiDiv_TEmetadata_all <- read.csv("widiv_TE_variation_matrix_revision_v1.csv") 
#oder protocol
LTRs_all = WiDiv_TEmetadata_all %>% filter(order == "LTR" & LTR_age != "NA")



#LTR age vs population frequency heatmap 
new_color = c("#8000FFFF","#00FFFFFF","#80FF00FF","#FF0000FF")
p <- ggplot(LTRs_all, aes(prop_present,LTR_age))
pdf("loess_regression_LTR_age_frequency.pdf",width =8, height = 5)
heatmap_pav = p + geom_bin2d(bins=150) +scale_fill_gradientn(limits = c(0,40),na.value = "red",labels = c("0","10","20","30",">=40"), breaks = seq(0,40,by=10), colours=new_color) + 
  #stat_smooth(aes(colour=("Regression"),method="gam",se = FALSE))+ theme_bw()
  #stat_smooth(method="loess", colour= "yellow", se = FALSE) + 
  labs(x = "Population Frequency", y = "LTR Similarity", fill='Count') + theme_classic() + theme(legend.position="top") +
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12))

# testing other color 

library(rasterVis)
library(httr)
library(viridis)
p <- ggplot(LTRs_all, aes(prop_present,LTR_age))
p + geom_bin2d(bins=150) + scale_fill_viridis(option = "plasma",limits = c(0,40),na.value = "yellow",labels = c("0","10","20","30","   40>=")) +
  stat_smooth(method="loess", colour= "yellow", se = FALSE) + 
  labs(x = "Population Frequency", y = "LTR Similarity", fill='Count') + theme_classic() + theme(legend.position="top") +
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12))



##########
# add age group and frequency group 
LTRs_all$class_stack<-ifelse(LTRs_all$prop_present <= 0.2,rr2<-"0-20%",
                                 ifelse(LTRs_all$prop_present>0.2 &LTRs_all$prop_present<=0.8,rr2<-"20%-80%",
                                        ifelse(LTRs_all$prop_present>0.8 & LTRs_all$prop_present<=1, rr2<-"80%-100%",
                                               rr2<-"")))

LTRs_all$age_stack<-ifelse(LTRs_all$LTR_age <= 95,rr2<-"Low Similarity",
                               ifelse(LTRs_all$LTR_age > 95 &LTRs_all$LTR_age<= 99,rr2<-"Moderate Similarity",
                                      ifelse(LTRs_all$LTR_age>99 & LTRs_all$LTR_age<=100, rr2<-"High Similarity",
                                             rr2<-"")))


# check the count of fixed very young TE 
fixed_very_young <- LTRs_all %>% filter(age_stack == "High Similarity" & LTRs_all$prop_present == 1)
dim(fixed_very_young)
#180 of very young TE are fixed 
write.csv(fixed_very_young,file="very_young_fixed_LTR.csv")

# count the number of fixed very young TE within 5kb of a gene
as.data.frame(table(fixed_very_young$genomic_loc))  %>% write.csv("fixed_very_young_genomic_location.csv")
# saved as the table fixed_very_young_TE_5kb_gene_frequency.csv with 33 that are within 5kb of the gene

# stats 

very_young_high_freq <- LTRs_all %>% filter(age_stack == "High Similarity" & class_stack =="80%-100%")

very_young_high_freq_count = as.data.frame(count(very_young_high_freq, 'family'))

# fam size of TEs 



WiDiv_TEmetadata_all_fam_size = WiDiv_TEmetadata_all %>% filter(LTR_age != "NA") %>% select(family)
df <- count(WiDiv_TEmetadata_all_fam_size,'family')

colnames(df) <- c("family","fam_size")
above_20 <- left_join(very_young_high_freq_count,df)



#
# LTR low similarity

LTRs_low_similarity = LTRs_all %>% filter(age_stack == "Low Similarity")
dim(LTRs_low_similarity)
table(LTRs_low_similarity$class_stack)
#   0-20%  20%-80% 80%-100% 
# 11760    45119    25461 

45119/82340 
#########split by TE age group 

#classify it into different TE age group 


table(LTRs_all$age_stack)
dim(LTRs_all)

summary_df <- LTRs_all %>%
  group_by(age_stack, class_stack) %>%
  summarise(n = n())

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

#break down into different frequency group then by LTR similarity 
LTR_low <- subset(LTRs_all, LTRs_all$class_stack == "0-20%")
LTRs_all$age_stack <- factor(LTRs_all$age_stack,levels =c("Low Similarity","Moderate Similarity","High Similarity"),ordered = TRUE)
nlabels_low <- table(LTR_low$age_stack)
LTR_low_violin <-ggplot(LTR_low, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_low),vjust = 3,hjust=-0.3) +
  labs( x = "Low Frequency LTR", y = "LTR Similarity (%)", colour='LTR Similarity Group') + theme_classic() +
  scale_color_manual(values=c("#1E88E5", "#FFC107", "#004D40")) + theme(axis.text.x=element_blank(),
                                                                        axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))

# LTR at moderate frequency

LTR_moderate <- subset(LTRs_all, LTRs_all$class_stack == "20%-80%")
nlabels_moderate<- table(LTR_moderate$age_stack)
LTR_moderate$age_stack <- factor(LTR_moderate$age_stack,levels =c("Low Similarity","Moderate Similarity","High Similarity"),ordered = TRUE)
LTR_moderate_violin <-ggplot(LTR_moderate, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal()  + stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_moderate),vjust = 3,hjust=-0.3) +
  labs( x = "Moderate Frequency LTR", y = "LTR Similarity (%)", colour='LTR Similarity Group') + theme_classic() + 
  scale_color_manual(values=c("#1E88E5", "#FFC107", "#004D40")) +   theme(axis.text.x=element_blank(),
                                                                          axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))


# LTR at high frequency
LTR_high <- subset(LTRs_all, LTRs_all$class_stack == "80%-100%")
LTR_high$age_stack <- factor(LTR_high$age_stack,levels =c("Low Similarity","Moderate Similarity","High Similarity"),ordered = TRUE)
nlabels_high<- table(LTR_high$age_stack)
LTR_high_violin <-ggplot(LTR_high, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() + stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_high),vjust = 3,hjust=-0.3) +
  labs( x = "High Frequency LTR", y = "LTR Similarity (%)", colour='LTR Similarity Group') + theme_classic() + 
  scale_color_manual(values=c("#1E88E5", "#FFC107", "#004D40")) + theme(axis.text.x=element_blank(),
                                                                        axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))


grid_arrange_shared_legend(LTR_low_violin,LTR_moderate_violin,LTR_high_violin, nrow = 3,ncol=1) 



# stats number for the manuscript LTR_high
fixed_high_frequency_fixed <- LTR_high %>% filter(prop_present =="1" & age_stack == "Very Young")
as.data.frame((table(fixed_high_frequency_fixed$genomic_loc)))

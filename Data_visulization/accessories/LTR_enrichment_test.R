setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
# LTR age 
WiDiv508_refB73_TEmetadata_age <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv") 
#oder protocol
B73_LTRs_all <- subset(WiDiv508_refB73_TEmetadata_age, order == "LTR")
dim(B73_LTRs_all)
B73_LTRs_age <- subset(B73_LTRs_all,B73_LTRs_all$LTR_age != "NA")
dim(B73_LTRs_age)
table(B73_LTRs_age$order)


#LTR age vs population frequency heatmap 
new_color = c("#8000FFFF","#00FFFFFF","#80FF00FF","#FF0000FF")
p <- ggplot(B73_LTRs_age, aes(prop_present,LTR_age))
p + geom_bin2d(bins=150) +scale_fill_gradientn(limits = c(0,40),na.value = "red",labels = c("0","10","20","30",">=40"), breaks = seq(0,40,by=10), colours=new_color) + 
  #stat_smooth(aes(colour=("Regression"),method="gam",se = FALSE))+ theme_bw()
  stat_smooth(method="loess", colour= "yellow", se = FALSE) + 
  labs(x = "Population Frequency", y = "LTR Similarity", fill='Count') + theme_classic() +
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12))


###############
# add age group and frequency group 
B73_LTRs_age$class_stack<-ifelse(B73_LTRs_age$prop_present <= 0.2,rr2<-"0-20%",
                             ifelse(B73_LTRs_age$prop_present>0.2 &B73_LTRs_age$prop_present<=0.8,rr2<-"20%-80%",
                                    ifelse(B73_LTRs_age$prop_present>0.8 & B73_LTRs_age$prop_present<=1, rr2<-"80%-100%",
                                           rr2<-"")))

B73_LTRs_age$age_stack<-ifelse(B73_LTRs_age$LTR_age <= 95,rr2<-"Old",
                           ifelse(B73_LTRs_age$LTR_age > 95 &B73_LTRs_age$LTR_age<= 99,rr2<-"Young",
                                  ifelse(B73_LTRs_age$LTR_age>99 & B73_LTRs_age$LTR_age<=100, rr2<-"Very Young",
                                         rr2<-"")))

##########
# check the count of fixed very young TE 
fixed_very_young <- subset(B73_LTRs_age,B73_LTRs_age$age_stack == "Very Young" & B73_LTRs_age$prop_present == 1)
dim(fixed_very_young)
#259 of very young TE are fixed 
write.csv(fixed_very_young,file="~/Desktop/TE paper/Fix_ambig/very_young_fixed_LTR.csv")

# count the number of fixed very young TE within 5kb of a gene
as.data.frame(table(fixed_very_young$genomic_loc)) 
# saved as the table fixed_very_young_TE_5kb_gene_frequency.csv with 33 that are within 5kb of the gene
 
#########
#split by TE age group 

#classify it into different TE age group 


table(B73_LTRs_age$age_stack)
dim(B73_LTRs_age)

summary_df <- B73_LTRs_age %>%
  group_by(age_stack, class_stack) %>%
  summarise(n = n())

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

#break down into different frequency group then by LTR similarity 
B73_LTR_low <- subset(B73_LTRs_age, B73_LTRs_age$class_stack == "0-20%")
B73_LTR_low$age_stack <- factor(B73_LTR_low$age_stack,levels =c("Old","Young","Very Young"),ordered = TRUE)

nlabels_low <- table(B73_LTR_low$age_stack)
B73_LTR_low_violin <-ggplot(B73_LTR_low, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_low),vjust = 3,hjust=-0.3) +
  labs( x = "Low Frequency LTR", y = "LTR Similarity (%)", colour='TE Age Group') + theme_classic() +
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF")) + theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
                                                                        axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))

# LTR at moderate frequency
B73_LTR_moderate <- subset(B73_LTRs_age, B73_LTRs_age$class_stack == "20%-80%")
B73_LTR_moderate$age_stack <- factor(B73_LTR_moderate$age_stack,levels =c("Old","Young","Very Young"),ordered = TRUE)
nlabels_moderate<- table(B73_LTR_moderate$age_stack)
B73_LTR_moderate_violin <-ggplot(B73_LTR_moderate, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal()  + stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_moderate),vjust = 3,hjust=-0.3) +
  labs( x = "Moderate Frequency LTR", y = "LTR Similarity (%)", colour='TE Age Group') + theme_classic() + 
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +   theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
                                                                          axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))


# LTR at high frequency
B73_LTR_high <- subset(B73_LTRs_age, B73_LTRs_age$class_stack == "80%-100%")

B73_LTR_high$age_stack <- factor(B73_LTR_high$age_stack,levels =c("Old","Young","Very Young"),ordered = TRUE)
nlabels_high<- table(B73_LTR_high$age_stack)
B73_LTR_high_violin <-ggplot(B73_LTR_high, aes(x=age_stack, y=LTR_age, color=age_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() + stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_high),vjust = 3,hjust=-0.3) +
  labs( x = "High Frequency LTR", y = "LTR Similarity (%)", colour='TE Age Group') + theme_classic() + 
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF")) + theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
                                                                        axis.text.y = element_text(color = "black", size = 12, face = "plain"), text = element_text(size=12))


grid_arrange_shared_legend(B73_LTR_low_violin,B73_LTR_moderate_violin,B73_LTR_high_violin, nrow = 3,ncol=1) 



###############################################Fisher Enrichment Test##########################

library(plyr)
# get set of high frequency LTR that are very young
Very_young_TE_high_frequency <- subset(B73_LTR_high,B73_LTR_high$age_stack == "Very Young")
dim(Very_young_TE_high_frequency) # 1854
#get a count of the very young TE at high frequency by family 
Very_young_TE_high_count <- count(Very_young_TE_high_frequency,"family") 
names(Very_young_TE_high_count)[2] <- "Very_Young_High_Frequency"
dim(Very_young_TE_high_count)
write.csv(Very_young_TE_high_count, file="Very_young_TE_fam_freq.csv")

#get total number for LTR at high frequency 
Very_young_TE_high_count["Total_very_young_LTR_high_freq"] = "1854"

head(Very_young_TE_high_count)
#all very young TE fam frequency at  level 
dim(B73_LTRs_age) # 177073
All_TE_All_Freq_count <- count(B73_LTRs_age,"family")
head(All_TE_All_Freq_count)

names(All_TE_All_Freq_count)[1]
names(All_TE_All_Freq_count)[2] <- "All_Frequency"
All_TE_All_Freq_count["total_LTR"] = "177073"

# Join matrix 
Join_young_TE_frequency <- left_join(All_TE_All_Freq_count, Very_young_TE_high_count)
na.omit(Join_young_TE_frequency)
write.csv(Join_young_TE_frequency,"Join_young_TE_frequency_20_above.csv")

# count 20 and above 

above_20 <- Join_young_TE_frequency %>% filter(Total_very_young_LTR_high_freq !=20) 

# fisher exact test 
greater_at_high_frequency <- apply(above_20, 1, 
                         function(x) {
                           tbl <- matrix(as.numeric(x[2:5]), ncol=2, byrow=T)
                           fisher.test(tbl, alternative="less")$p.value
                         })
names(greater_at_high_frequency)[1] <- "greater_at_high_frequency"
greater_at_high_frequency <- as.data.frame(greater_at_high_frequency)


final_enrichment_list <- cbind(above_20,greater_at_high_frequency)
write.csv(final_enrichment_list, file="~/Desktop/TE paper/Revision_MS/Table S3. enrichment_very_young_high_freq_all_list.csv")


# total enriched (494 familes, no size restriction)
table(final_enrichment_list$greater_at_high_frequency < 0.01)
#FALSE  TRUE 
#463    31 


# fam size 20 or above
#total number of TE familes 
fam_size_20<- subset(final_enrichment_list,final_enrichment_list$All_Frequency >= 20)
dim(fam_size_20)
write.csv(fam_size_20, file="~/Desktop/TE paper/Fix_ambig/enrichment_fam_size_20.csv")

fam_size <- subset(final_enrichment_list,final_enrichment_list$All_Frequency)
dim(fam_size)
write.csv(fam_size, file="~/Desktop/TE paper/Fix_ambig/enrichment_fam_size.csv")


#enriched at p_value 0.01 
fam_size_20_0.01 <- subset(final_enrichment_list,final_enrichment_list$All_Frequency >= 20 & final_enrichment_list$greater_at_high_frequency < 0.01)
dim(fam_size_20_0.01)
write.csv(fam_size_20_0.01, file="enrichment_fam_size_20_0.01.csv")


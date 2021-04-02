library(reshape2)
library(Rmisc)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(dplyr)
library(scales)
library(lemon)
library(grid)
library(gridExtra)

setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")

full_matrix <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv")
LTR_matrix <- full_matrix %>% filter(order == "LTR" & LTR_age != "NA")
arm_info <- read.csv("B73_LTR_chr_region_info.txt",sep = ',',header=FALSE)
colnames(arm_info) <- c("TE","chr_region")

join_B73_LTR_matrix <- left_join(arm_info,LTR_matrix)

join_B73_LTR_matrix$class_stack<-ifelse(join_B73_LTR_matrix$prop_present <= 0.2,rr2<-"0-20%",
                                 ifelse(join_B73_LTR_matrix$prop_present>0.2 &join_B73_LTR_matrix$prop_present<=0.8,rr2<-"20%-80%",
                                        ifelse(join_B73_LTR_matrix$prop_present>0.8 & join_B73_LTR_matrix$prop_present<=1, rr2<-"80%-100%",
                                               rr2<-"")))

join_B73_LTR_matrix$age_stack<-ifelse(join_B73_LTR_matrix$LTR_age <= 95,rr2<-"Old",
                               ifelse(join_B73_LTR_matrix$LTR_age > 95 &join_B73_LTR_matrix$LTR_age<= 99,rr2<-"Young",
                                      ifelse(join_B73_LTR_matrix$LTR_age>99 & join_B73_LTR_matrix$LTR_age<=100, rr2<-"Very Young",
                                             rr2<-"")))

#(1)
#fixed_high_frequency_very_young_LTR
very_young_high_freq = join_B73_LTR_matrix %>% filter(age_stack =="Very Young" & class_stack =="80%-100%"& prop_present =="1")
table(very_young_high_freq$chr_region)
#arm        pericentromeric 
#26                 142  

#(2)
# get TE family informaiton 
fixed_element_family <- count(very_young_high_freq,"family")
# total 168 families

# subset dataframe only include TE familes that have fixed_element above
subset_all_TE_fam <- join_B73_LTR_matrix[is.element(join_B73_LTR_matrix$family, very_young_high_freq$family),]
fixed_element_family_all <- count(subset_all_TE_fam,"family")
table(subset_all_TE_fam$chr_region)

#arm        pericentromeric 
#35742             30836 

# so for the rest TE elements in the TE families we are looking, they are not fixed. 
#arm        pericentromeric 
#35742-26 = 35716          30836-142= 30694


#(3)
# genome-wide fixed LTR
fixed_LTR_genome_wide <- join_B73_LTR_matrix %>% filter(prop_present =="1")
table(fixed_LTR_genome_wide$chr_region)
#arm        pericentromeric 
#1927              6230

#(4) genome-wide that are very young
very_young_genome_wide <- join_B73_LTR_matrix %>% filter(age_stack == "Very Young")
table(very_young_genome_wide$chr_region)

#arm          pericentromeric 
#9212              6445

#(5)
table(join_B73_LTR_matrix$chr_region)
#arm          pericentromeric 
#57473             51495


write.csv(join_B73_LTR_matrix,file="~/Desktop/TE_chromsome_distribution/LTR_age_chr_region.csv")

# fisher test 
fixed_high_freq_vs_all = as.data.frame(rbind(c("arm","Fixed_Very_Young",26),
                                             c("pericentromeric","Fixed_Very_Young",142),
                                             c("arm","Unfixed_very_young_fam",35716),
                                             c("pericentromeric","Unfixed_very_young_fam",30694),
                                             c("arm","All_fixed_LTR",1927),
                                             c("pericentromeric","All_fixed_LTR",6230),
                                             c("arm","Genome-wide all_very_young",9212),
                                             c("pericentromeric","Genome-wide all_very_young",6445),
                                             c("arm","All_LTRs",57473),
                                             c("pericentromeric","All_LTRs",51495)))
                                          
colnames(fixed_high_freq_vs_all) <- c("region","type","count")


write.csv(fixed_high_freq_vs_all,file = "chr_region_LTR_age_stack_plot.csv")

stack <- read.csv("chr_region_LTR_age_stack_plot.csv")

stack$type <- factor(stack$type, levels = c("Fixed_Very_Young","All_fixed_LTR","Unfixed_very_young_fam","Genome-wide all_very_young","All_LTRs"))


ggplot(data=stack,aes(x=type,y=count,fill=region)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) + 
  theme(text = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=12, angle=0,vjust = 0.5),
        axis.text.y = element_text(color="black", size=10, angle=30),axis.title.y = element_text(size = 12)) + guides(fill=guide_legend(title="Recombination Frequency"))+
        ylab("Percentage") + xlab("") + scale_fill_manual(values=c("#5F4B8BFF", "#E69A8DFF"),labels=c("High Recombination","Low Recombination")) +
  scale_x_discrete(labels=c("Fixed_Very_Young" = "Fixed High Similarity LTRs", "All_fixed_LTR" = "All fixed LTRs", "Unfixed_very_young_fam"="Unfixed High Similarity LTRs",
                            "Genome-wide all_very_young" = "Genome-wide All High Similarity LTRs","All_LTRs"="All LTRs")) + coord_flip() + theme(legend.position="center")  



  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) #+
  scale_x_discrete(breaks = c(3,6,9))


fisher.test(fixed_high_freq_vs_all,alternative = "greater")
#p-value = 5.437e-05 


unfixed_vs_all <- rbind(c(30694,35716),
                     c(51495,57473))

fisher.test(unfixed_vs_all,alternative = "less")

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

all_set <- read.csv("~/Desktop/TE paper/File_S5_final.txt", sep = "\t")
arm_info <- read.csv("~/Desktop/TE_chromsome_distribution/B73_LTR_chr_region_info.txt")
LTR_matrix <- subset(all_set,all_set$ambig_cat == "low_ambig" & all_set$LTR_age !="NA")

join_B73_LTR_matrix <- left_join(LTR_matrix,arm_info)

join_B73_LTR_matrix$class_stack<-ifelse(join_B73_LTR_matrix$prop_present <= 0.2,rr2<-"0-20%",
                                 ifelse(join_B73_LTR_matrix$prop_present>0.2 &join_B73_LTR_matrix$prop_present<=0.8,rr2<-"20%-80%",
                                        ifelse(join_B73_LTR_matrix$prop_present>0.8 & join_B73_LTR_matrix$prop_present<=1, rr2<-"80%-100%",
                                               rr2<-"")))

join_B73_LTR_matrix$age_stack<-ifelse(join_B73_LTR_matrix$LTR_age <= 95,rr2<-"Old",
                               ifelse(join_B73_LTR_matrix$LTR_age > 95 &join_B73_LTR_matrix$LTR_age<= 99,rr2<-"Young",
                                      ifelse(join_B73_LTR_matrix$LTR_age>99 & join_B73_LTR_matrix$LTR_age<=100, rr2<-"Very Young",
                                             rr2<-"")))


working_set <- na.omit(join_B73_LTR_matrix)
#(1)
#fixed_high_frequency_very_young_LTR
very_young_high_freq <- subset(working_set,working_set$age_stack =="Very Young" & working_set$class_stack =="80%-100%"& working_set$prop_present =="1" )
table(very_young_high_freq$chr_region)
#arm        pericentromeric 
#58               184  

#(2)
# get TE family informaiton 
fixed_element_family <- count(very_young_high_freq,"family")
# total 117 families 
# subset dataframe only include TE familes that have fixed_element above
subset_all_TE_fam <- join_B73_LTR_matrix[is.element(join_B73_LTR_matrix$family, fixed_element_family$family),]
fixed_element_family_all <- count(subset_all_TE_fam,"family")
table(subset_all_TE_fam$chr_region)

#arm        pericentromeric 
#32050             27973 

# so for the rest TE elements in the TE families we are looking
#arm        pericentromeric 
#32050-58 = 31992          27973-184= 27789


#(3)
# genome-wide fixed LTR
fixed_LTR_genome_wide <- subset(working_set,working_set$prop_present =="1")
table(fixed_LTR_genome_wide$chr_region)
#arm        pericentromeric 
#3309              8858

#(4) genome-wide that are young
very_young_genome_wide <- subset(working_set,working_set$age_stack == "Very Young")
table(very_young_genome_wide$chr_region)

#arm          pericentromeric 
#8404             5842

#(5)
table(working_set$chr_region)
#arm          pericentromeric 
#48837             45617


write.csv(join_B73_LTR_matrix,file="~/Desktop/TE_chromsome_distribution/LTR_age_chr_region.csv")

# fisher test 
fixed_high_freq_vs_all = as.data.frame(rbind(c("arm","Fixed_Very_Young",58),
                                             c("pericentromeric","Fixed_Very_Young",184),
                                             c("arm","Unfixed_very_young_fam",31992),
                                             c("pericentromeric","Unfixed_very_young_fam",27789),
                                             c("arm","All_fixed_TE",3309),
                                             c("pericentromeric","All_fixed_TE",8858),
                                             c("arm","Genome-wide all_very_young",8404),
                                             c("pericentromeric","Genome-wide all_very_young",5842),
                                             c("arm","All_TEs",48837),
                                             c("pericentromeric","All_TEs",45617)))
                                          
colnames(fixed_high_freq_vs_all) <- c("region","type","count")


write.csv(fixed_high_freq_vs_all,file = "~/Desktop/TE_chromsome_distribution/stack_plot.csv")

stack <- read.csv("~/Desktop/TE_chromsome_distribution/stack_plot.csv")

stack$type <- factor(stack$type, levels = c("Fixed_Very_Young","All_fixed_TE","Unfixed_very_young_fam","Genome-wide all_very_young","All_TEs"))


ggplot(data=stack,aes(x=type,y=count,fill=region)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) + 
  theme(text = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=12, angle=80,vjust = 0.5),
        axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12)) + guides(fill=guide_legend(title="Recombination Frequency"))+
        ylab("Percentage") + xlab("") + scale_fill_manual(values=c("#5F4B8BFF", "#E69A8DFF"),labels=c("High Recombination","Low Recombination")) + 
  scale_x_discrete(breaks=c("All LTR","Genome-wide Young LTR", "Genome-wide Fixed LTR","Fixed Very Young LTR","Unfixed LTR in Fixed Very Young LTR Family"))
  



  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) #+
  scale_x_discrete(breaks = c(3,6,9))


fisher.test(fixed_high_freq_vs_all,alternative = "greater")
#p-value = 5.437e-05 


unfixed_vs_all <- rbind(c(27789,31992),
                     c(45617,48837))

fisher.test(unfixed_vs_all,alternative = "less")

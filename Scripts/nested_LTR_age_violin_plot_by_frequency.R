library(ggplot2)
library(ggthemes)
library(extrafont)
library(dplyr)
library(scales)
library(lemon)
library(grid)
library(gridExtra)
library(tidyverse)

setwd("~/Desktop/TE paper/bioRxiv_figure_nested_vs_outer_TE_age_frequency/")

WiDiv508_refB73_nest_out<- read.csv("WiDiv508_B73_nested.v.outer_present.txt",sep = "\t")

table(WiDiv508_refB73_nest_out$nestedTE_ambigcat)
table(WiDiv508_refB73_nest_out$outerTE_ambigcat)

# filter outer TE with presence more than 0.95 
outer_TE_0.95 <- subset(WiDiv508_refB73_nest_out, WiDiv508_refB73_nest_out$outerTE_prop_present > 0.95)
dim(outer_TE_0.95)

outer_TE_0.95$outerTE_name_for_order <- substr(outer_TE_0.95$outerTE_name, 0, 3)
table(outer_TE_0.95$outerTE_name_for_order)
outer_TE_0.95$nested_TE_name_for_order <- substr(outer_TE_0.95$nestedTE_name, 0, 3)


outer_TE_0.95$nested_TE_order<-ifelse(outer_TE_0.95$nested_TE_name_for_order == "DHH",rr2 <-"Helitron",
                                     ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTA",rr2 <-"TIR",
                                            ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTC", rr2<-"TIR",
                                                   ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTH", rr2<-"TIR", 
                                                          ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTM", rr2<-"TIR", 
                                                                 ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTT", rr2<-"TIR",
                                                                        ifelse(outer_TE_0.95$nested_TE_name_for_order == "DTX", rr2<-"TIR",
                                                                               ifelse(outer_TE_0.95$nested_TE_name_for_order == "RIL", rr2<-"LINE",
                                                                                      ifelse(outer_TE_0.95$nested_TE_name_for_order == "RIT", rr2<-"LINE",
                                                                                          ifelse(outer_TE_0.95$nested_TE_name_for_order == "RLC", rr2<-"LTR",
                                                                                              ifelse(outer_TE_0.95$nested_TE_name_for_order == "RLG", rr2<-"LTR",
                                                                                                  ifelse(outer_TE_0.95$nested_TE_name_for_order == "RLX", rr2<-"LTR",
                                                                                                         ifelse(outer_TE_0.95$nested_TE_name_for_order == "RST", rr2<-"SINE",
                                                                                                      rr2<-"")))))))))))))


table(outer_TE_0.95$nested_TE_order)

outer_TE_0.95$outer_TE_order <- ifelse(outer_TE_0.95$outerTE_name_for_order == "DHH",rr2 <-"Helitron",
                                       ifelse(outer_TE_0.95$outerTE_name_for_order == "DTA",rr2 <-"TIR",
                                              ifelse(outer_TE_0.95$outerTE_name_for_order == "DTC", rr2<-"TIR",
                                                     ifelse(outer_TE_0.95$outerTE_name_for_order == "DTH", rr2<-"TIR", 
                                                            ifelse(outer_TE_0.95$outerTE_name_for_order == "DTM", rr2<-"TIR", 
                                                                   ifelse(outer_TE_0.95$outerTE_name_for_order == "DTT", rr2<-"TIR",
                                                                          ifelse(outer_TE_0.95$outerTE_name_for_order == "DTX", rr2<-"TIR",
                                                                                 ifelse(outer_TE_0.95$outerTE_name_for_order == "RIL", rr2<-"LINE",
                                                                                        ifelse(outer_TE_0.95$outerTE_name_for_order == "RIT", rr2<-"LINE",
                                                                                               ifelse(outer_TE_0.95$outerTE_name_for_order == "RLC", rr2<-"LTR",
                                                                                                      ifelse(outer_TE_0.95$outerTE_name_for_order == "RLG", rr2<-"LTR",
                                                                                                             ifelse(outer_TE_0.95$outerTE_name_for_order == "RLX", rr2<-"LTR",
                                                                                                                    ifelse(outer_TE_0.95$outerTE_name_for_order == "RST", rr2<-"SINE",
                                                                                                                           rr2<-"")))))))))))))
table(outer_TE_0.95$outer_TE_order)
# getting dataset with nested TE from LTR 

nested_LTR_0.95 <- subset(outer_TE_0.95,outer_TE_0.95$nested_TE_order == "LTR")
dim(nested_LTR_0.95)

# adding LTR age data via left join using nested_LTR_0.95 nestedTE_name as the query 
age_matrix <- read.csv("age_matrix_for_left_join.csv")

matrix_with_age <- na.omit(left_join(nested_LTR_0.95[-1,],age_matrix))
dim(matrix_with_age)

# define frequency group using nestedTE_prop_present 
matrix_with_age$nested_class_stack<-ifelse(matrix_with_age$nestedTE_prop_present <= 0.2,rr2<-"0-20%",
                            ifelse(matrix_with_age$nestedTE_prop_present>0.2 &matrix_with_age$nestedTE_prop_present<=0.8,rr2<-"20%-80%",
                                   ifelse(matrix_with_age$nestedTE_prop_present>0.8 & matrix_with_age$nestedTE_prop_present<=1, rr2<-"80%-100%",
                                          rr2<-"")))

table(matrix_with_age$nested_class_stack)

write.csv(matrix_with_age,file = "matrix_with_age.csv")
matrix_with_age$nested_class_stack <- as.factor(matrix_with_age$nested_class_stack)

# this gives out the number for each sample size 
as.data.frame(matrix_with_age$nested_class_stack )
summary(as.data.frame(matrix_with_age$nested_class_stack ))
#664
#5055 
#12090
nlabels_all_outer_TE <- table(matrix_with_age$nested_class_stack)
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}


# 3 nested LTR has outer element LINEs, so the total number will not match with the number when split them by class
All_outer_TE_order <-ggplot(matrix_with_age, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_all_outer_TE),vjust = 5,hjust=-0.4,size = 5) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
          axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        text = element_text(size=12))
                                                                                                                           




# all_outer_TE with outer TE LTR/TIR/Helitron

major_outer <- subset(matrix_with_age,matrix_with_age$outer_TE_order != "LINE" & matrix_with_age$outer_TE_order != "SINE")  
dim(major_outer)
dim(matrix_with_age)
nlabels_all_outer_TE <- table(major_outer$nested_class_stack)

table(major_outer$outer_TE_order)

All_major_outer_TE_order <-ggplot(major_outer, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_all_outer_TE),vjust = 5,hjust=-0.4,size = 5) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        text = element_text(size=12))

# split outer TE by the order 

# LTR 
outer_LTR <- subset(matrix_with_age,matrix_with_age$outer_TE_order == "LTR")

nlabels_outer_LTR <- table(outer_LTR$nested_class_stack)
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# make sure the outer elements 
LTR_outer_TE_order <-ggplot(outer_LTR, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_outer_LTR),vjust = 12,hjust=-0.1) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + theme(text = element_text(size = 12)) + 
theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
      axis.text.y = element_text(color = "black", size = 12, face = "plain"),
      text = element_text(size=12))

#TIR 
outer_TIR <- subset(matrix_with_age,matrix_with_age$outer_TE_order == "TIR")
table(outer_TIR$outer_TE_order)
nlabels_outer_TIR <- table(outer_TIR$nested_class_stack)
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

TIR_outer_TE_order <-ggplot(outer_TIR, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_outer_TIR),vjust = 5,hjust=-0.6) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + theme(text = element_text(size = 12)) + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        text = element_text(size=12))

#Helitron
outer_Helitron <- subset(matrix_with_age,matrix_with_age$outer_TE_order == "Helitron")
nlabels_outer_Helitron <- table(outer_Helitron$nested_class_stack)
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

Helitron_outer_TE_order <-ggplot(outer_Helitron, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_outer_Helitron),vjust = 5,hjust=-0.4) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + theme(text = element_text(size = 12)) + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        text = element_text(size=12))



grid_arrange_shared_legend(LTR_outer_TE_order, TIR_outer_TE_order, Helitron_outer_TE_order, nrow = 1,ncol=3) 

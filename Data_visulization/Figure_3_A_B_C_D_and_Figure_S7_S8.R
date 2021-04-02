setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")

TE_meta_metrix <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv")

#B73 nested TE inner outer information
B73_inner_outer_base <- read.csv(file = "Figure3/B73.structuralTEv2.2018.12.20.filteredTE_nested.and.outerTEs.txt",sep = "\t",header = TRUE)
colnames(B73_inner_outer_base) <- c("nestedTE_name","outerTE_name")
#getting prop_present for each TE
TE_freq_nested <- TE_meta_metrix %>% select("TE","prop_present","order")
colnames(TE_freq_nested) <- c("nestedTE_name","nestedTE_prop_present","nested_order")


TE_freq_outer <- TE_meta_metrix %>% select("TE","prop_present","order")
colnames(TE_freq_outer) <- c("outerTE_name","outerTE_prop_present","outer_order")

# make join matrix 

nested_outer_matrix <- left_join(B73_inner_outer_base,TE_freq_outer,by="outerTE_name") %>% left_join(TE_freq_nested,by="nestedTE_name")

# remove rows that contains NA 
nested_outer_matrix <- na.omit(nested_outer_matrix) 
# 168444 TEs in this analysis 






a_correlation_nested = nested_outer_matrix %>% ggplot(aes(nestedTE_prop_present,outerTE_prop_present)) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("Proportion of Genotypes the Nested TE is Presented in") + ylab("Proportion of Genotypes the Outer TE is Presented in") + 
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())


# stats for manuscript 
nested_outer_matrix

# panel b
nested_outer_matrix$delta_freq <- nested_outer_matrix$outerTE_prop_present - nested_outer_matrix$nestedTE_prop_present
below005 <- nested_outer_matrix %>% filter(delta_freq <0 &delta_freq >= -0.05 )
above000 <- nested_outer_matrix %>% filter(delta_freq >=0)
dim(above000)
dim(nested_outer_matrix)
#120372 
# total
# 168444

b_histgram_frequency <- nested_outer_matrix %>% ggplot(aes(delta_freq)) + geom_histogram(color="black", fill="grey") + theme_classic() + 
  xlab("Outer element frequency – Nested element frequency") + ylab("Frequency") + 
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())
b_histgram_frequency


##add in similarity 

outer_LTR_similarity <- read.csv(file = "~/Desktop/TE paper/Revision/joining_matrix/LTR_similarity_YQ.txt",sep = "\t",header=FALSE)
colnames(outer_LTR_similarity) <- c("outerTE_name","outer_age")
nestedTE_name_LTR_similarity <- read.csv(file = "~/Desktop/TE paper/Revision/joining_matrix/LTR_similarity_YQ.txt",sep = "\t",header=FALSE)
colnames(nestedTE_name_LTR_similarity) <- c("nestedTE_name","nestedTE_name_age")

# total 42534 TE are used for the analysis 
# delta frequency has been calculated in previous section, only need to do age difference here

# Figure 3C 

# for this analysis, it doesn't matter what outer TEs are, as long as it has a nested LTR

outer_095 <- nested_outer_matrix %>% filter(outerTE_prop_present >0.95)
outer_095_nested_LTR_age <- left_join(outer_095,nestedTE_name_LTR_similarity,by="nestedTE_name")
outer_095_nested_LTR_age <- na.omit(outer_095_nested_LTR_age) 
outer_095_nested_LTR_age$nested_class_stack<-ifelse(outer_095_nested_LTR_age$nestedTE_prop_present <= 0.2,rr2<-"0-20%",
                                           ifelse(outer_095_nested_LTR_age$nestedTE_prop_present>0.2 &outer_095_nested_LTR_age$nestedTE_prop_present<=0.8,rr2<-"20%-80%",
                                                  ifelse(outer_095_nested_LTR_age$nestedTE_prop_present>0.8 & outer_095_nested_LTR_age$nestedTE_prop_present<=1, rr2<-"80%-100%",
                                                         rr2<-"")))



outer_095_nested_LTR_age %>% ggplot(aes(nested_class_stack,nestedTE_name_age)) + geom_violin(trim=TRUE) + geom_boxplot(width=0.1) +theme_classic() +
  xlab("Nested Element Population Frequency") + ylab("LTR Similarity (%)") + theme(axis.title.y = element_text(size=12),
                                                                                   axis.title.x = element_text(size=12),
                                                                                   axis.text.y = element_text(size=12), 
                                                                                   axis.text.x = element_text(size=12),
                                                                                   legend.text = element_text(size=12),
                                                                                   legend.title=element_blank()) 


table(outer_095_nested_LTR_age$nested_class_stack)
# 3 nested LTR has outer element LINEs, so the total number will not match with the number when split them by class
All_outer_TE_order <-ggplot(matrix_with_age, aes(x=nested_class_stack, y=LTR_age, color=nested_class_stack)) +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1) + theme_minimal() +#stat_summary(fun.data = give.n, geom = "text",label = paste("n =", nlabels_all_outer_TE),vjust = 5,hjust=-0.4,size = 5) + 
  labs( x = "Nested Element Population Frequency", y = "LTR Similarity (%)", colour='TE Frequency') + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),
        text = element_text(size=12))


# figure 3D 
LTR_age_freq <- left_join(LTR_only,outer_LTR_similarity)  %>% left_join(nestedTE_name_LTR_similarity,by="nestedTE_name")
LTR_age_freq <- na.omit(LTR_age_freq) 
LTR_age_freq$delta_age <- LTR_age_freq$outer_age - LTR_age_freq$nestedTE_name_age

# stats for manuscript 
similarity_inner_high <- LTR_age_freq %>% filter(delta_age <0)
# 35780 

# total 
42534
Figure_3D = LTR_age_freq %>% ggplot(aes(delta_freq,delta_age)) +  theme_classic() + 
  geom_rect(aes(xmin=0, xmax=1, ymin=0,ymax=-15), alpha=0.01, fill="yellow") + 
  geom_rect(aes(xmin=-1, xmax=1, ymin=0,ymax=15), alpha=0.01, fill="grey") + 
  geom_rect(aes(xmin=-1, xmax=0, ymin=0,ymax=-15), alpha=0.01, fill="grey") + 
  geom_point(alpha = 1/10) + 
  geom_vline(xintercept = 0, color = "white", size=0.5) + 
  geom_hline(yintercept = 0, color = "white", size=0.5) + 
  xlab("Frequency:Outer Element-Nested Element") + ylab("LTR Similarity: Outer Element-Nested Element") + 
  theme_classic()
  

Figure_3D + theme(axis.title.y = element_text(size=12),
                    axis.title.x = element_text(size=12),
                    axis.text.y = element_text(size=12), 
                    axis.text.x = element_text(size=12),
                    legend.text = element_text(size=12),
                    legend.title=element_blank()) 


# Figure S7

outer_nested_age_all <- left_join(nested_outer_matrix,outer_LTR_similarity)  %>% left_join(nestedTE_name_LTR_similarity,by="nestedTE_name")

# outer element > 0.95 
outer_fixed_all = outer_nested_age_all %>% filter(outerTE_prop_present > 0.95)
nested_age = outer_fixed_all%>% filter(nestedTE_name_age != "NA")

nested_age$nested_class_stack<-ifelse(nested_age$nestedTE_prop_present <= 0.2,rr2<-"0-20%",
                                                    ifelse(nested_age$nestedTE_prop_present>0.2 &nested_age$nestedTE_prop_present<=0.8,rr2<-"20%-80%",
                                                           ifelse(nested_age$nestedTE_prop_present>0.8 & nested_age$nestedTE_prop_present<=1, rr2<-"80%-100%",
                                                                  rr2<-"")))



no_line$outer_order = nested_age %>% filter(outer_order != "LINE" ) %>%  filter(nested_class_stack != "NA" )

no_line$outer_order <- factor(no_line$outer_order,levels=c("LTR","TIR","Helitron"))

ggplot(no_line,aes(x = nested_class_stack, y=nestedTE_name_age))+ geom_violin(trim=TRUE) + geom_boxplot(width=0.1) +
  facet_wrap(no_line$outer_order,scales = "free",strip.position = "bottom") +theme_classic()   + 
  scale_x_discrete(name="Nested TE Population Frequency") + scale_y_continuous("LTR Similarity (%)")+theme_classic() + 
  theme(strip.background =element_rect(fill="grey"),text = element_text(size = 12)) 

check = no_line %>% filter(outer_order == "TIR")
table(check$nested_class_stack)

check = no_line %>% filter(outer_order == "LTR")
table(check$nested_class_stack)

check = no_line %>% filter(outer_order == "Helitron")
table(check$nested_class_stack)



# Figure S8 
age_compare <- outer_nested_age_all %>% filter(outer_age != "NA" & nestedTE_name_age != "NA") %>% drop_na()
age_compare %>% 
  ggplot(aes(nestedTE_name_age,outer_age))  +
  geom_point(alpha = 1/10) + theme_classic() + 
  xlim(85, 100) + ylim(85, 100) + 
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank()) + xlab("Nested Element LTR Similarity (%)") + ylab("Outer Element LTR Similarity (%)") + 
  geom_abline(intercept = 0, slope = 1,colour='red')

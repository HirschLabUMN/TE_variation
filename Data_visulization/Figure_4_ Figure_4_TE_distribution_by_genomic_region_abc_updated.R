library(reshape2)
library(Rmisc)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(dplyr)
library(scales)
library(lemon)
library(grid)
library(lattice)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# read final matrix 
df <- read.csv("~/Desktop/Widiv_TE_variation_matrix_revision_v1_fmt.csv")
# features to loop over and subsetting each time: 
feature_loop <- c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")

# subset matrix with LTRs that has age information 
# starting with Helitron 
Helitron_subset <- df %>% filter(order=="Helitron")

# for loop function
ploting_matrix=NULL
for (i in feature_loop){
  TE_feature <- Helitron_subset %>% filter(genomic_loc== i)
  TE_feature$group <- "Matrix"
  # this information will be used to feed in the random sampling 
  TE_number_count = nrow(TE_feature)
  # random sampling the total number of the TE for each feature and repeat 100 times
  drawing_mean=NULL
  for (j in 1:100) {
    random_drawing = sample_n(Helitron_subset, TE_number_count)
    feature_info <- as.data.frame(random_drawing[1,])
    feature_info$prop_present <- mean(random_drawing$prop_present) #replace the prop_present value using the mean of the drawing
    drawing_mean <- rbind(drawing_mean,feature_info) # adding each mean value from the 100 iternation to the matrix 
  }
  # adding a group class to the matrix 
  #extract line info to replace value for visulization. This prop_present is now standing for the mean frequency drawing from this random selection + bootstrapping 
  drawing_mean$group <- "random_sampling"
  # replace the genomic loc to the same feature used in the loop so can be plot together
  drawing_mean$genomic_loc <- paste0(i)
  feature_subset <- rbind(TE_feature,drawing_mean)
  ploting_matrix <- rbind(ploting_matrix,feature_subset)
}

whole_set <- ploting_matrix  
whole_set$genomic_class <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                 ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                        ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                               ifelse(whole_set$genomic_loc == "intron", rr2<-"within_gene", 
                                                      ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                             ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                    ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                           ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"downstream",
                                                                                  ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"upstream",
                                                                                         ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"downstream",
                                                                                                ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"upstream",
                                                                                                       ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                              ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                                     ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"upstream",
                                                                                                                            rr2<-""))))))))))))))
whole_set$genomic_class_fine <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                      ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                             ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                                    ifelse(whole_set$genomic_loc == "intron", rr2<-"intron", 
                                                           ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                                  ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                         ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                                ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"near_gene",
                                                                                       ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"near_gene",
                                                                                              ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"near_gene",
                                                                                                     ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"near_gene",
                                                                                                            ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                   ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                          ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"near_gene",
                                                                                                                                 rr2<-""))))))))))))))



whole_set$genomic_loc <- factor(whole_set$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
whole_set$genomic_class <- factor(whole_set$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 
whole_set$group <- factor(whole_set$group, levels = c("Matrix","random_sampling")) 


#rename this set as Helitron_wholeset to store data for plotting. 
Helitron_wholeset <- whole_set 


############ moving on to LTR set 
LTR_subset <- df %>% filter(order=="LTR") %>% filter(LTR_age != "NA")
table(LTR_subset$genomic_loc)
# for loop function
ploting_matrix=NULL
for (i in feature_loop){
  TE_feature <- LTR_subset %>% filter(genomic_loc== i)
  TE_feature$group <- "Matrix"
  # this information will be used to feed in the random sampling 
  TE_number_count = nrow(TE_feature)
  # random sampling the total number of the TE for each feature and repeat 100 times
  drawing_mean=NULL
  for (j in 1:100) {
    random_drawing = sample_n(LTR_subset, TE_number_count)
    feature_info <- as.data.frame(random_drawing[1,])
    feature_info$prop_present <- mean(random_drawing$prop_present) #replace the prop_present value using the mean of the drawing
    drawing_mean <- rbind(drawing_mean,feature_info) # adding each mean value from the 100 iternation to the matrix 
  }
  # adding a group class to the matrix 
  #extract line info to replace value for visulization. This prop_present is now standing for the mean frequency drawing from this random selection + bootstrapping 
  drawing_mean$group <- "random_sampling"
  # replace the genomic loc to the same feature used in the loop so can be plot together
  drawing_mean$genomic_loc <- paste0(i)
  feature_subset <- rbind(TE_feature,drawing_mean)
  ploting_matrix <- rbind(ploting_matrix,feature_subset)
}

whole_set <- ploting_matrix  
whole_set$genomic_class <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                 ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                        ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                               ifelse(whole_set$genomic_loc == "intron", rr2<-"within_gene", 
                                                      ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                             ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                    ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                           ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"downstream",
                                                                                  ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"upstream",
                                                                                         ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"downstream",
                                                                                                ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"upstream",
                                                                                                       ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                              ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                                     ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"upstream",
                                                                                                                            rr2<-""))))))))))))))
whole_set$genomic_class_fine <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                      ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                             ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                                    ifelse(whole_set$genomic_loc == "intron", rr2<-"intron", 
                                                           ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                                  ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                         ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                                ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"near_gene",
                                                                                       ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"near_gene",
                                                                                              ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"near_gene",
                                                                                                     ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"near_gene",
                                                                                                            ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                   ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                          ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"near_gene",
                                                                                                                                 rr2<-""))))))))))))))



whole_set$genomic_loc <- factor(whole_set$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
whole_set$genomic_class <- factor(whole_set$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 
whole_set$group <- factor(whole_set$group, levels = c("Matrix","random_sampling")) 
#rename this set as Helitron_wholeset to store data for plotting. 
LTR_wholeset <- whole_set 


########
# TIR section 
TIR_subset <- df %>% filter(order=="TIR")

# for loop function
ploting_matrix=NULL
for (i in feature_loop){
  TE_feature <- TIR_subset %>% filter(genomic_loc== i)
  TE_feature$group <- "Matrix"
  # this information will be used to feed in the random sampling 
  TE_number_count = nrow(TE_feature)
  # random sampling the total number of the TE for each feature and repeat 100 times
  drawing_mean=NULL
  for (j in 1:100) {
    random_drawing = sample_n(TIR_subset, TE_number_count)
    feature_info <- as.data.frame(random_drawing[1,])
    feature_info$prop_present <- mean(random_drawing$prop_present) #replace the prop_present value using the mean of the drawing
    drawing_mean <- rbind(drawing_mean,feature_info) # adding each mean value from the 100 iternation to the matrix 
  }
  # adding a group class to the matrix 
  #extract line info to replace value for visulization. This prop_present is now standing for the mean frequency drawing from this random selection + bootstrapping 
  drawing_mean$group <- "random_sampling"
  # replace the genomic loc to the same feature used in the loop so can be plot together
  drawing_mean$genomic_loc <- paste0(i)
  feature_subset <- rbind(TE_feature,drawing_mean)
  ploting_matrix <- rbind(ploting_matrix,feature_subset)
}

whole_set <- ploting_matrix  
whole_set$genomic_class <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                 ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                        ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                               ifelse(whole_set$genomic_loc == "intron", rr2<-"within_gene", 
                                                      ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                             ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                    ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                           ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"downstream",
                                                                                  ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"upstream",
                                                                                         ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"downstream",
                                                                                                ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"upstream",
                                                                                                       ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                              ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"downstream",
                                                                                                                     ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"upstream",
                                                                                                                            rr2<-""))))))))))))))
whole_set$genomic_class_fine <-ifelse(whole_set$genomic_loc == "3prime_UTR",rr2 <-"within_gene",
                                      ifelse(whole_set$genomic_loc == "5prime_UTR",rr2 <-"within_gene",
                                             ifelse(whole_set$genomic_loc == "exon", rr2<-"within_gene",
                                                    ifelse(whole_set$genomic_loc == "intron", rr2<-"intron", 
                                                           ifelse(whole_set$genomic_loc == "TE_encompassed_by_gene", rr2<-"within_gene", 
                                                                  ifelse(whole_set$genomic_loc == "TE_encompassing_gene", rr2<-"TE_encompassing_gene",
                                                                         ifelse(whole_set$genomic_loc == "TE_intergenic", rr2<-"TE_intergenic",
                                                                                ifelse(whole_set$genomic_loc == "TE_1.5kb_downstream", rr2<-"near_gene",
                                                                                       ifelse(whole_set$genomic_loc == "TE_10.5kb_upstream", rr2<-"near_gene",
                                                                                              ifelse(whole_set$genomic_loc == "TE_1kb_downstream", rr2<-"near_gene",
                                                                                                     ifelse(whole_set$genomic_loc == "TE_1kb_upstream", rr2<-"near_gene",
                                                                                                            ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                   ifelse(whole_set$genomic_loc == "TE_5.10kb_downstream", rr2<-"near_gene",
                                                                                                                          ifelse(whole_set$genomic_loc == "TE_5.1kb_upstream", rr2<-"near_gene",
                                                                                                                                 rr2<-""))))))))))))))



whole_set$genomic_loc <- factor(whole_set$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
whole_set$genomic_class <- factor(whole_set$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 
whole_set$group <- factor(whole_set$group, levels = c("Matrix","random_sampling")) 


#rename this set as Helitron_wholeset to store data for plotting. 
TIR_wholeset <- whole_set 


########
# all dataset are LTR_wholeset, TIR_wholeset, and Helitron_wholeset


# Figure visulization
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
scaleFUN <- function(x) sprintf("%.2f", x)
color_levels = c("#B18CD9","#FFB85F","#FBCCD1","#00AAA0")


# Helitron
p_Helitron<- ggplot(Helitron_wholeset, aes(x=genomic_loc, y=prop_present,fill=genomic_class,alpha=group)) + scale_alpha_discrete(range = c(1, 0.2)) +
  geom_boxplot(width=0.4) + #facet_wrap(~genomic_class,scales = "free") +
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("Helitron Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_Helitron_join <- p_Helitron + scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("upstream","within_gene","TE_encompassing_gene","downstream"))     


# TIR 
p_TIR<- ggplot(TIR_wholeset, aes(x=genomic_loc, y=prop_present,fill=genomic_class,alpha=group)) + scale_alpha_discrete(range = c(1, 0.2)) +
  geom_boxplot(width=0.4) + #facet_wrap(~genomic_class,scales = "free") +
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("TIR Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_TIR_join <- p_TIR + scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("upstream","within_gene","TE_encompassing_gene","downstream"))     

# LTR 
p_LTR <- ggplot(LTR_wholeset, aes(x=genomic_loc, y=prop_present,fill=genomic_class,alpha=group)) + scale_alpha_discrete(range = c(1, 0.2)) +
  geom_boxplot(width=0.4) + #stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")+
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("LTR Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_LTR_join <- p_LTR + scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("upstream","within_gene","TE_encompassing_gene","downstream"))     

ggarrange(p_Helitron_join,p_TIR_join,p_LTR_join,ncol=1)

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

setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
whole_set <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv")
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

LTR_subset <- subset(whole_set,whole_set$order == "LTR")
dim(LTR_subset)
# [1] 190373     16

# remove TEs from intergenic regions and with age information
working_LTR_age <- subset(LTR_subset,LTR_subset$genomic_loc != "TE_intergenic" & LTR_subset$LTR_age !="NA")
dim(working_LTR_age)
# [1] 37838    16

working_LTR_age$age_stack<-ifelse(working_LTR_age$LTR_age > 99,rr2<-"High Similarity",
                                  ifelse(working_LTR_age$LTR_age>95 & working_LTR_age$LTR_age<=99,rr2<-"Moderate Similarity",
                                         ifelse(working_LTR_age$LTR_age <=95, rr2<-"Low Similarity",
                                                rr2<-"")))

dfwc_age <- summarySEwithin(working_LTR_age, measurevar="prop_present", withinvars="genomic_loc",
                            idvar="order", na.rm=FALSE, conf.interval=.95)


as.data.frame(dfwc_age)

dfwc_age$genomic_loc <- factor(dfwc_age$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
scaleFUN <- function(x) sprintf("%.2f", x)


working_LTR_age$genomic_loc <- factor(working_LTR_age$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
working_LTR_age$genomic_class <- factor(working_LTR_age$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 
table(working_LTR_age$genomic_class)
library(RColorBrewer)
# LTR box plot 
color_levels = c("#B18CD9","#FFB85F","#FBCCD1","#00AAA0")
p_LTR <- ggplot(working_LTR_age, aes(x=genomic_loc, y=prop_present,fill=genomic_class)) +
  geom_boxplot(width=0.4) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("LTR Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_LTR +   scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("Upstream","Gene Body","TE Encompassing Gene","Downstream"))     
p1 <- p_LTR +   scale_fill_manual()     



# stack bar plot (Figure 4D)
working_LTR_age$genomic_loc <- factor(working_LTR_age$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream"))
working_LTR_age$age_stack <- factor(working_LTR_age$age_stack, levels = c("Low Similarity","Moderate Similarity","High Similarity"))

LTR_age_stack_pct <- ggplot(working_LTR_age, aes(factor(genomic_loc),fill = age_stack)) +
  geom_bar(stat="count", position = "fill")  



LTR_age_stack_pct2 <- LTR_age_stack_pct + theme(axis.text.x = element_text(color="black", 
                                                                           size=12, angle=80,vjust = 0.5),
                                                axis.text.y = element_text(color="black", 
                                                                           size=12, angle=0)) + 
  theme(legend.position="top") + ylab("LTR Similarity Group Composition") + labs(fill="LTR Similarity Group") + xlab("") +
  scale_fill_manual(values=c("#1E88E5", "#FFC107", "#004D40")) + 
  theme(text = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(color="black", size=12, angle=80,vjust = 0.5),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12)) +
  scale_x_discrete(breaks=c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream"),labels=c("5-10 kb upstream","1-5 kb upstream","0-1 kb upstream","5' UTR","exon","intron","3' UTR","TE encompassed by gene","TE encompassing gene", "0-1 kb downstream","1-5 kb downstream","5-10 kbdownstream"))


########enrichment test
###enchriment of very young TE in the gene region 
#total number of very young lTR 
table(working_LTR_age$order)
#total number of LTR with age info 
#[1] 30407    19
# new 37838 
table(working_LTR_age$genomic_class)
#near_gene TE_encompa sing_gene          within_gene 
#27766                 1380                 1261  

# new  downstream TE_encompassing_gene             upstream          within_gene 
# 15687                 1779                18732                 1640 

#LTR in the gene region and summarize the age 
LTR_in_gene <- subset(working_LTR_age,working_LTR_age$genomic_class == "within_gene")
table(LTR_in_gene$age_stack)
#Old      Young Very Young 
#464        499        298 

# new 
# Low Similarity Moderate Similarity     High Similarity 
# 628                 635                 377 
#LTR in any region
table(working_LTR_age$age_stack)
#Old      Young Very Young 
#14135      11690       4582 


# new      Low Similarity Moderate Similarity     High Similarity 
# 18402               14239                5197 


## very young TE in the gene region (including intron 298), total 1261
## very young TE in any gene region 4582, total 30407
# enrichment test
LTR_age_enrichment <- rbind(
  c(377,1640),
  c(5197,37838)
)

fisher.test(LTR_age_enrichment,alternative="greater")


###for genemic_class  TE_encompassing_gene - a total of 1881
LTR_encompassing_gene<- subset(working_LTR_age,working_LTR_age$genomic_class == "TE_encompassing_gene")
table(working_LTR_age$genomic_class)
#near_gene TE_encompassing_gene          within_gene 
#27766                 1380                 1261 
# new downstream TE_encompassing_gene             upstream          within_gene 
# 15687                 1779                18732                 1640 
table(LTR_encompassing_gene$age_stack)
#Old      Young Very Young 
#1103        673        105 

# new Low Similarity Moderate Similarity     High Similarity 
# 1057                 622                 100 

LTR_age_encompassing_gene <- rbind(
  c(100,1779),
  c(5197,37838)
)

fisher.test(LTR_age_encompassing_gene,alternative="less")





# helitron
#Helitron 
Helitron_subset0 <- subset(whole_set,whole_set$order == "Helitron")
Helitron_subset <- subset(Helitron_subset0,Helitron_subset0$genomic_loc != "TE_intergenic")
# total number of Helitron
table(Helitron_subset$order)
#Helitron     LINE      LTR     SINE      TIR 
#6171        0        0        0        0 
# helitron that are in within
table(Helitron_subset$genomic_class)
#near_gene TE_encompassing_gene          within_gene 
#4356                 1540                  275 


dfwc_Helitron_subset <- summarySEwithin(Helitron_subset, measurevar="prop_present", withinvars="genomic_loc",
                                        idvar="order", na.rm=FALSE, conf.interval=.95)


dfwc_Helitron_subset$genomic_loc <- factor(dfwc_Helitron_subset$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream"))
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
scaleFUN <- function(x) sprintf("%.2f", x)
Helitron_frequency_ave <- ggplot(dfwc_Helitron_subset, aes(x=genomic_loc, y=prop_present, group=1)) +
  #geom_line() +
  geom_errorbar(width=.1, aes(ymin=prop_present-se, ymax=prop_present+se))  + 
  geom_text(aes(y = 0.5,label = paste0("\nN=",N))) + ylab("Helitron Population Frequency") +
  geom_point(shape=21, size=2, fill="white")   + scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.10),labels = fmt_dcimals(2)) + 
  theme(text = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=16, angle=0),axis.title.y = element_text(size = 16))


table(Helitron_subset$genomic_class)


table(Helitron_subset$genomic_class_fine)
#intron            near_gene TE_encompassing_gene          within_gene 
#130                 4356                 1540                  145 

table(Helitron_subset$order)


#violin_helitron
Helitron_subset$genomic_loc <- factor(Helitron_subset$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
Helitron_subset$genomic_class <- factor(Helitron_subset$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 


# hilitron box plot 
color_levels = c("#B18CD9","#FFB85F","#FBCCD1","#00AAA0")
#color_levels = c("white","white","white","white")
p_Helitron <- ggplot(Helitron_subset, aes(x=genomic_loc, y=prop_present,fill=genomic_class)) +
  geom_boxplot(width=0.4) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("Helitron Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_Helitron +   scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("Upstream","Gene Body","TE Encompassing Gene","Downstream"))     



# TIR 
TIR_subset0 <- subset(whole_set,whole_set$order == "TIR")
TIR_subset <- subset(TIR_subset0,TIR_subset0$genomic_loc != "TE_intergenic")

# total number of TIR
table(TIR_subset$order)
#Helitron     LINE      LTR     SINE      TIR 
#0        0        0        0    42598 
# helitron that are in TIR
table(TIR_subset$genomic_class)
#           near_gene TE_encompassing_gene          within_gene 
#                   73                 4874 
dfwc_TIR_subset <- summarySEwithin(TIR_subset, measurevar="prop_present", withinvars="genomic_loc",
                                   idvar="order", na.rm=FALSE, conf.interval=.95)

dfwc_TIR_subset$genomic_loc <- factor(dfwc_TIR_subset$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream"))


fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
scaleFUN <- function(x) sprintf("%.2f", x)
TIR_frequency_ave <- ggplot(dfwc_TIR_subset, aes(x=genomic_loc, y=prop_present, group=1)) +
  #geom_line() +
  geom_errorbar(width=.1, aes(ymin=prop_present-se, ymax=prop_present+se))  + 
  geom_text(aes(y = 0.5,label = paste0("\nN=",N))) + ylab("TIR Population Frequency") +
  geom_point(shape=21, size=2, fill="white")  + scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + xlab("") + 
  theme(text = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(color="black", 
                                                                                                                  size=16, angle=80,vjust = 0.5),axis.text.y = element_text(color="black", size=16, angle=0),axis.title.y = element_text(size = 16)) +
  scale_x_discrete(breaks=c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream"),labels=c("5-10 kb upstream","1-5 kb upstream","0-1 kb upstream","5' UTR","exon","intron","3' UTR","TE encompassed by gene","TE encompassing gene", "0-1 kb downstream","1-5 kb downstream","5-10 kbdownstream")) 

#check intron TE number 
table(TIR_subset$genomic_class_fine)
table(TIR_subset$genomic_loc)

#intron            near_gene TE_encompassing_gene          within_gene 
#3151                37651                   73                 1723 
table(TIR_subset$order)

#Helitron     LINE      LTR     SINE      TIR 
#0        0        0        0    42598 



#TIR box plot 
TIR_subset$genomic_loc <- factor(TIR_subset$genomic_loc, levels = c("TE_10.5kb_upstream","TE_5.1kb_upstream","TE_1kb_upstream","5prime_UTR","exon","intron","3prime_UTR","TE_encompassed_by_gene","TE_encompassing_gene", "TE_1kb_downstream","TE_1.5kb_downstream","TE_5.10kb_downstream")) 
TIR_subset$genomic_class <- factor(TIR_subset$genomic_class, levels = c("upstream","within_gene","TE_encompassing_gene","downstream")) 
table(TIR_subset$genomic_class)
library(RColorBrewer)

color_levels = c("#B18CD9","#FFB85F","#FBCCD1","#00AAA0")
#color_levels = c("white","white","white","white")
p_TIR <- ggplot(TIR_subset, aes(x=genomic_loc, y=prop_present,fill=genomic_class)) +
  geom_boxplot(width=0.4) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1.0, by = 0.20),labels = fmt_dcimals(2)) + ylab("TIR Population Frequency") + 
  theme(text = element_text(size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y = element_text(color="black", size=12, angle=0),axis.title.y = element_text(size = 12),legend.position="none")  
p_TIR +   scale_fill_manual(values=c(color_levels[1],color_levels[2],color_levels[3],color_levels[4]),name = "",labels = c("Upstream","Gene Body","TE Encompassing Gene","Downstream"))     



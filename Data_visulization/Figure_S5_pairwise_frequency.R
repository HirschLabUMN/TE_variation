setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
library(tidyverse)
library(ggpubr)
library(gridExtra)

freq_4ref <- read.csv("~/Desktop/TE paper/Revision/joining_matrix/freq_matrix_four_refs.csv")
table(freq_4ref$pop_freq_B73)

# frequency comparision 
## B73 vs Mo17
B73_Mo17 = freq_4ref %>% select("pop_freq_B73","pop_freq_Mo17") %>% drop_na() %>% 
  ggplot(aes(pop_freq_B73,pop_freq_Mo17))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (B73 Reference Genome)") + ylab("TE Population Frequency (Mo17 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

## B73 vs PH207
B73_PH207 = freq_4ref %>% select("pop_freq_B73","pop_freq_PH207") %>% drop_na() %>% 
  ggplot(aes(pop_freq_B73,pop_freq_PH207))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (B73 Reference Genome)") + ylab("TE Population Frequency (PH207 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

## B73 vs W22
B73_W22 = freq_4ref %>% select("pop_freq_B73","pop_freq_W22") %>% drop_na() %>% 
  ggplot(aes(pop_freq_B73,pop_freq_W22))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (B73 Reference Genome)") + ylab("TE Population Frequency (W22 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

## Mo17 vs W22
Mo17_W22 = freq_4ref %>% select("pop_freq_Mo17","pop_freq_W22") %>% drop_na() %>% 
  ggplot(aes(pop_freq_Mo17,pop_freq_W22))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (Mo17 Reference Genome)") + ylab("TE Population Frequency (W22 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

## Mo17 vs PH207
Mo17_PH207 = freq_4ref %>% select("pop_freq_Mo17","pop_freq_PH207") %>% drop_na() %>% 
  ggplot(aes(pop_freq_Mo17,pop_freq_PH207))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (Mo17 Reference Genome)") + ylab("TE Population Frequency (PH207 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

## PH207 vs W22
PH207_W22 = freq_4ref %>% select("pop_freq_PH207","pop_freq_W22") %>% drop_na() %>% 
  ggplot(aes(pop_freq_PH207,pop_freq_W22))  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.6, label.y = 0.05) + geom_point(alpha = 1/50) + theme_classic() + 
  xlab("TE Population Frequency (PH207 Reference Genome)") + ylab("TE Population Frequency (W22 Reference Genome)") + 
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank())

grid.arrange(B73_Mo17,B73_Mo17,B73_W22,Mo17_PH207,Mo17_W22,PH207_W22)

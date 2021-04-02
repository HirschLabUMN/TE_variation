setwd("/Users/yinjieqiu/Desktop/TE paper/Revision/mean_genome_coverage/")

library(ggridges)

refB73_509 <- read.csv(file = "refB73v4_mean_cov_allbams.txt_509.txt",sep = '\t',header=FALSE)
refMo17_509 <- read.csv(file = "refMo17_mean_cov_allbams.txt_509.txt",sep = '\t',header=FALSE)
refPH207_509 <- read.csv(file = "refPH207_mean_cov_allbams.txt_509.txt",sep = '\t',header=FALSE)
refW22_509 <- read.csv(file = "refW22v12_mean_cov_allbams.txt_509.txt",sep = '\t',header=FALSE)


coverage_matrix <- rbind(refB73_509,refMo17_509,refPH207_509,refW22_509)
coverage_matrix %>% 
  ggplot(aes(x=V3,y=V2)) + geom_density_ridges()  +
  xlab("Realized Genome-Wide Mean Coverage (x)") + ylab("Refence Genomes") + theme_classic()+ 
  theme(axis.text = element_text(size = 12),
        axis.title =element_text(size = 12)) + 
  geom_vline(xintercept = 15, linetype="dotted",color = "blue", size=1) + 
  geom_vline(xintercept = 30, linetype="dotted",color = "blue", size=1)

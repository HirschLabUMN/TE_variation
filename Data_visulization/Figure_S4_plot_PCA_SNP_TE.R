library(lemon)
library(tidyverse)

# TE based PCA 
setwd("~/Desktop/TE paper/Revision_MS/data_visulization/Figure_S4/")
widiv_info <- read.csv("~/Desktop/TE paper/Revision_MS/data_visulization/widiv_line_info.csv")
colnames(widiv_info) <- c("taxa", "Maize Type")


# PCA TE 
pca_TE <- read_table2("plink_pca_TE.eigenvec", col_names = FALSE)
eigenval_TE <- scan("plink_pca_TE.eigenval")
pve_TE = eigenval_TE/sum(eigenval_TE)*100
pve_TE <- data.frame(PC = 1:10, pve_TE = eigenval_TE[1:10]/sum(eigenval_TE[1:10])*100)
pca_TE <- pca_TE[,-1]
names(pca_TE)[1] <- "taxa"
names(pca_TE)[2:ncol(pca_TE)] <- paste0("PC", 1:(ncol(pca_TE)-1))
data_TE <- left_join(pca_TE,widiv_info)
str(data_TE)
summary(data_TE$`Maize Type`)
unique(data_TE$`Maize Type`)
print(levels(data_TE$`Maize Type`))
TE_PCA = ggplot(data_TE) +
  geom_point(aes(x = PC1, y = PC2, color = `Maize Type`), alpha = 0.5, size = 2) +
  theme_classic() +
  xlab(paste0("PC1 (", signif(pve_TE$pve_TE[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_TE$pve_TE[2], 3), "%)")) + ggtitle("Widiv TE PCA") +
  theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank()) + theme(legend.position = "none")



# SNP full PCA 
pca <- read_table2("plink_pca_no_filter.eigenvec", col_names = FALSE)
eigenval <- scan("plink_pca_no_filter.eigenval")


pve <- data.frame(PC = 1:10, pve = eigenval[1:10]/sum(eigenval[1:10])*100)
pca <- pca[,-1]
names(pca)[1] <- "taxa"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
data_SNP <- left_join(pca,widiv_info)
str(data_SNP)
summary(data_SNP$`Maize Type`)
unique(data_SNP$`Maize Type`)
print(levels(data_SNP$`Maize Type`))
#write.csv(data_SNP,"data_SNP_no_filter.csv")

#data_SNP <- read.csv("data_SNP_no_filter.csv")
SNP_PCA_full= ggplot(data_SNP) +
  geom_point(aes(x = PC1, y = PC2, color = `Maize Type`), alpha = 0.5, size = 2) +
  theme_classic() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle("Widiv SNP PCA") +
  theme(axis.title.y = element_text(size=10), axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title=element_blank()) +theme(legend.position = "none")


grid_arrange_shared_legend(SNP_PCA_full,TE_PCA,nrow=1)


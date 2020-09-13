#analyze B73 gene expression data to identify representative tissues 
library(tidyverse)

#Why am I  doing this? 
#Have gene/TE expression data for many tissues and developemntal timepoints. 
#A lot of this information will be redundant
#I want to find representative tissues 

############################################################################
#graph TE expression data
#B73 TEs
TE_rpm_tissue <- read_tsv("~/Documents/TE_project/gene_te_expression_data/element_RPM_B73_development_11Oct19.txt")
TE_rpm_tissue$superfam <- substr(TE_rpm_tissue$TE_ID, 1,3)
TE_rpm_pca <- prcomp(TE_rpm_tissue[,c(2:364)], center = T, scale = T)

fviz_pca_ind(TE_rpm_pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = FALSE)
fviz_pca_var(TE_rpm_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), geom = c("point","text"), repeal = F)

#read in file with variable contirbutions to PC1 and 2
element_rpm_pca_var <- read_tsv("~/Documents/TE_project/gene_te_expression_data/element_rpm_pca_var.txt")
ggplot(element_rpm_pca_var, aes(x = pc_dim1, y = pc_dim2, color = partial_tissue_name)) + geom_point()

#try to plot TEs colored by superfamily 
TE_rpm_pca_ind <- get_pca_ind(TE_rpm_pca)
TE_rpm_pca_ind_coord <- TE_rpm_pca_ind$coord
plot(TE_rpm_pca_ind$coord[,1], TE_rpm_pca_ind$coord[,2], xlim = c(0,500), ylim = c(-150,50), main = "TE expression PCA")
#add TE superfam name
#this assumes original order of TEs in file was preserved, which I think it is
TE_rpm_pca_ind_coord$superfam <- TE_rpm_tissue$superfam

TE_pca_indiv_dim12 <- cbind(TE_rpm_pca_ind_coord[,1], TE_rpm_pca_ind_coord[,2], TE_rpm_tissue$superfam)
ggplot(TE_rpm_pca_ind_coord, aes(x = TE_rpm_pca_ind_coord[,1], y = TE_rpm_pca_ind_coord[,2], color = TE_rpm_pca_ind_coord[,300])) + geom_point()
ggplot(TE_pca_indiv_dim12, aes(x = as.numeric([,1]), y = as.numeric([,2]))) + geom_point()

############################################################################
#analyze gene expression data
library(factoextra)
#above is for graphing 

#B73_geneexp <- read_tsv("~/Documents/TE_project/gene_RPM_B73_development_20Nov19.txt")
#use prcomp to run pca on gene expression data
#goal: identify the tissues (columns) that contribute most to gene expression variation and remove redundant tissues

B73_geneexp_prcomp <- prcomp(B73_geneexp[,2:364], scale = T)
fviz_eig(B73_geneexp_prcomp)
#PC1 explains ~42% of the variation
#graph individuals/genes: 
fviz_pca_ind(B73_geneexp_prcomp, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = F, geom = "point")
#at a gross visual level the variation looks pretty similar to TE variation (a big clump + a few outilers on both axes)

#graph variables (tissue, this is what I'm more interested in)
fviz_pca_var(B73_geneexp_prcomp, col.var = "x", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = F, geom = "point")
#colored by contribution to the pc

fviz_pca_biplot(B73_geneexp_prcomp, col.var = "#2E9FDF", geom.var = "point", col.ind = "#696969")

#what variables/tissues are most important? 
#contribution of variables to dim-1/PC1
fviz_contrib(B73_geneexp_prcomp, choice = "var", axes = 1)
fviz_contrib(B73_geneexp_prcomp, choice = "var", axes = 2)
fviz_contrib(B73_geneexp_prcomp, choice = "var", axes = 1:3)

fviz_contrib(B73_geneexp_prcomp, choice = "var", axes = 1:10, top = 230)


B73_geneexp_pca_var <- get_pca_var(B73_geneexp_prcomp)

B73_geneexp_eigenval <- get_eigenvalue(B73_geneexp_prcomp)
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

#read in transposed gene expression data file
B73_geneexp_tp <- read_tsv("~/Documents/TE_project/gene_te_expression_data/gene_RPM_B73_development_20Nov19_col1tissue_names.txt")
#get mean gene expression values
B73_geneexp_tp_mean <- aggregate(B73_geneexp_tp[,4:36226], list(B73_geneexp_tp$tissue_name_full), mean)
#transpose mean gene expression data frame
#make first column row names 
rownames(B73_geneexp_tp_mean) <- B73_geneexp_tp_mean$Group.1
B73_geneexp_tp_mean$Group.1 <- NULL
#convert to matrix
B73_geneexp_mean_retranspose <- as.data.frame(t(as.matrix(B73_geneexp_tp_mean)))
#create tree
corr <- cor(na.omit(B73_geneexp_mean_retranspose))
dist <- as.dist(1-corr)
clust <- hclust(dist)
plot(clust)
#ggdendrogram(clust, rotate = F, size = 2)

library(stringr)
B73_tissue_name_split <- (str_split_fixed(B73_geneexp_tp_mean$Group.1, "_", 2))
#add column with just tissue name
B73_geneexp_tp_mean$tissue_name_split <- B73_tissue_name_split[,1]

#do pca on tissues
B73_mean_genexp_pca <- prcomp(B73_geneexp_tp_mean[,2:36224])
fviz_pca_ind(B73_mean_genexp_pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = F, geom = "point")
#create new dataframe 
B73_tissue_PC1PC2 <- data.frame(B73_mean_genexp_pca$x[,1:2])
B73_tissue_PC1PC2$tissue_name_full <- B73_geneexp_tp_mean$Group.1
B73_tissue_PC1PC2$tissue_type <- B73_geneexp_tp_mean$tissue_name_split
ggplot(B73_tissue_PC1PC2, aes(x = PC2, y = PC1, color = tissue_type)) + geom_point() + ggtitle("Tissue mean gene expression PCA")
ggsave(filename = "~/Documents/TE_project/gene_te_expression_data/B73_gene_tissue_exp_PCA.png", device = "png", width = 10, height = 6, units = "in")

#k-means clustering 
#rows are observations and columns are variables 
#scale RPM values 
B73_geneexp_mean_retranspose <- scale(B73_geneexp_mean_retranspose)
B73_mean_geneexp_kmeans <- kmeans(x = B73_geneexp_mean_retranspose, centers = 10, nstart = 25)
fviz_cluster(B73_mean_geneexp_kmeans, data = B73_geneexp_mean_retranspose)
B73_mean_geneexp_kmeans$centers
#I think above is what I would use
#get tissue that is closest to each center 
#what are those tissues? 
head(sort(B73_mean_geneexp_kmeans$centers[1,]))
#lots of overlap in tissues closest to each center 
#what does this mean? 


#distance <- get_dist(B73_geneexp_mean_retranspose)

#I think kmeans clustering needs to be done with the tissues as rows 
#nope, rows are observstions and columns are variables
B73_mean_geneexp_tissuerow_kmeans <- kmeans(x = B73_geneexp_tp_mean, centers = 10, nstart = 25)
fviz_cluster(B73_mean_geneexp_tissuerow_kmeans, data = B73_geneexp_tp_mean)

######################################################################################
#select tissues and look at average expression
B73_geneexp_mean_retranspose <- as.data.frame(t(as.matrix(B73_geneexp_tp_mean)))
#row names to column names 
library(data.table)
B73_geneexp_mean_retranspose <- setDT(B73_geneexp_mean_retranspose, keep.rownames = T)[]
B73_geneexp_select_tissues <- subset(B73_geneexp_mean_retranspose, select = c("rn","endosperm_16dap_DS", "anther_R1_DS", "leaf_18dap_DS", "leaf_sheath.V12_DZ", "root_primary.GH.6das_DS", "cob_immature.V18_DS", "tassel_meiotic.V18_DS", "internode_18dap_DS", "seed_12dap_DS", "embryo_16dap_DS"))
#B73_geneexp_select_tissues <- data.frame(B73_geneexp_select_tissues)
B73_geneexp_select_tissues$exp_count <- 0
B73_geneexp_select_tissues$exp_prop <- 0

#add count of #tissues with mean RPM > 1
for (i in 1:nrow(B73_geneexp_select_tissues)) {
	count_1 <- length(which(B73_geneexp_select_tissues[i,2:11] > 1))
	B73_geneexp_select_tissues[i,12] <- count_1
	prop_1 <- count_1/10
	B73_geneexp_select_tissues[i,13] <- prop_1
}
ggplot(B73_geneexp_select_tissues, aes(x = exp_prop)) + geom_histogram(bins = 10) + ggtitle("Proportion of 10 representative tissues a gene is expressed in\nexpressed = mean RPM > 1") + xlab("Proportion of tissues expressed")
ggsave(filename = "~/Documents/TE_project/gene_te_expression_data/B73_gene_tissue_exp_proportion.png", device = "png", width = 8, height = 5, units = "in")
write.table(B73_geneexp_select_tissues, file = "~/Documents/TE_project/gene_te_expression_data/B73_gene_select_tissue_exp.txt", row.names = F, sep = "\t")
##############################################
#now read in TE expression data (by family) and get mean expression and distribution
TE_fam_rpm_tissue <- read_tsv("~/Documents/TE_project/gene_te_expression_data/TE_RPM_matrix_REDO_11Oct19.txt")
#transpose to be able to calculate mean expression 
rownames(TE_fam_rpm_tissue) <- TE_fam_rpm_tissue$TE_famname
TE_fam_rpm_tissue$TE_famname <- NULL
#convert to matrix
TE_fam_rpm_transpose <- as.data.frame(t(as.matrix(TE_fam_rpm_tissue)))
TE_fam_rpm_transpose$tissue_name_full <- rownames(TE_fam_rpm_transpose)
#split based on tissue name 
library(stringr)
TE_fam_tissue_name_split <- (str_split_fixed(TE_fam_rpm_transpose$tissue_name_full, "_", 3))
TE_fam_tissue_name_split_rep <- (str_split_fixed(TE_fam_rpm_transpose$tissue_name_full, "_B", 2))
#add column with just tissue name
TE_fam_rpm_transpose$tissue_name_no_rep <- TE_fam_tissue_name_split_rep[,1]
TE_fam_rpm_transpose$tissue_name_split <- TE_fam_tissue_name_split[,1]

TE_fam_rpm_mean <- aggregate(TE_fam_rpm_transpose[,1:23129], list(TE_fam_rpm_transpose$tissue_name_no_rep), mean)
#re-transpose data frame to be able to select 10 tissues 
rownames(TE_fam_rpm_mean) <- TE_fam_rpm_mean$Group.1
TE_fam_rpm_mean$Group.1 <- NULL
#convert to matrix
TE_fam_rpm_mean_retranspose <- as.data.frame(t(as.matrix(TE_fam_rpm_mean)))
TE_fam_rpm_mean_retranspose$TE_fam_name <- rownames(TE_fam_rpm_mean_retranspose)
#select 10 tissues and get #present etc 
TEfam_exp_select_tissues <- subset(TE_fam_rpm_mean_retranspose, select = c("endosperm_16dap_DS", "anther_R1_DS", "leaf_18dap_DS", "leaf_sheath.V12_DZ", "root_primary.GH.6das_DS", "cob_immature.V18_DS", "tassel_meiotic.V18_DS", "internode_18dap_DS", "seed_12dap_DS", "embryo_16dap_DS"))
TEfam_exp_select_tissues$exp_count <- 0
TEfam_exp_select_tissues$exp_prop <- 0
TEfam_exp_select_tissues$te_fam_name <- rownames(TEfam_exp_select_tissues)
TEfam_exp_select_tissues$superfam <- substr(TEfam_exp_select_tissues$te_fam_name,1,3)

#add count of #tissues with mean RPM > 1
for (i in 1:nrow(TEfam_exp_select_tissues)) {
	count_1 <- length(which(TEfam_exp_select_tissues[i,1:10] > 1))
	TEfam_exp_select_tissues[i,11] <- count_1
	prop_1 <- count_1/10
	TEfam_exp_select_tissues[i,12] <- prop_1
}
ggplot(TEfam_exp_select_tissues, aes(x = exp_prop)) + geom_histogram(bins = 10) + ggtitle("Proportion of 10 representative tissues a TE family is expressed in\nexpressed = mean RPM > 1") + xlab("Proportion of tissues expressed") 
ggsave(filename = "~/Documents/TE_project/gene_te_expression_data/TEfam_tissue_exp_proportion.png", device = "png", width = 8, height = 5, units = "in")
write.table(TEfam_exp_select_tissues, file = "~/Documents/TE_project/gene_te_expression_data/TEfam_select_tissue_exp.txt", row.names = F, sep = "\t")


#end of script

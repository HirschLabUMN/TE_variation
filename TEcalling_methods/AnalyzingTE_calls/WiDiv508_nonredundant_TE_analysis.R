#TE project: figures and modeling for non-redundant TE dataset 

library(tidyverse)
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(plotrix)
library(data.table)

#read in non redundant data sets 
non_redundant_allTEs <- read_tsv(file = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/non_redundant_TEs_WiDiv508_meancounts_TE.metadata.txt")
non_redundant_TEs_LTRs <- read_tsv(file = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/non_redundant_TEs_WiDiv508_meancounts_TE.metadata_LTRs.txt")
non_redundant_TEs_within10kb <- read_tsv(file = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/non_redundant_TEs_WiDiv508_meancounts_TE.metadata_bp.togenes.upneg.txt")
non_redundant_TEs_rawcounts <- read_tsv(file = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/non_redundant_TEs_WiDiv508_counts_non0.txt")

#B73 counts vs Mo17 counts 
non_redundant_TEs_rawcounts_B73.Mo17 <- subset(non_redundant_TEs_rawcounts, B73.status == "present" & Mo17.status == "present" & B73.n_ambiguous < 128 & Mo17.n_ambiguous < 128)
ggplot(non_redundant_TEs_rawcounts_B73.Mo17, aes(x = B73.n_present, y = Mo17.n_present)) + geom_point(alpha = 0.02) + ggtitle("Number of genotypes a TE was called as 'Present' in between homologous TEs in B73 and Mo17")
ggplot(non_redundant_TEs_rawcounts_B73.Mo17, aes(x = as.numeric(B73.n_present/(B73.n_absent+B73.n_present)), y = as.numeric(Mo17.n_present/(Mo17.n_absent+Mo17.n_present)))) + geom_point(alpha = 0.02) + ggtitle("Proportion of genotypes a TE was called as 'Present' in between homologous TEs in B73 and Mo17 \n low ambiguous TEs only") + xlab("Proportion present in B73 reference genome TE") + ylab("Proportion present in Mo17 reference genome TE")
ggsave(filename = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/WiDiv508_proppresent_B73_v_Mo17.png", device = "png", width = 7.5, height = 8, units = "in")
ggsave(filename = "~/Dropbox/HirschLab_MaizeWork/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/WiDiv508_proppresent_B73_v_Mo17.svg", device = "svg", width = 7.5, height = 8, units = "in")

#process non-redundant dataset 
#choose only low ambig TEs and for TIRs, helitrons and LTRs 
non_redundant_TEmetadata_allTEs_low.ambig <- subset(non_redundant_allTEs, ambig_cat == "low_ambig" & order == "Helitron" | ambig_cat == "low_ambig" & order == "TIR" | ambig_cat == "low_ambig" & order == "LTR")
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig <- subset(non_redundant_TEs_within10kb, ambig_cat == "low_ambig" & order == "Helitron" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "TIR" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "LTR" & distance_to_gene != "te_part_ingene")
non_redundant_TEs_LTRs_low.ambig <- subset(non_redundant_TEs_LTRs, ambig_cat == "low_ambig")

#get mean proportion present by genomic location bins of all TEs
non_redundant_TEmetadata_low.ambig_sum_genomic_loc <- non_redundant_TEmetadata_allTEs_low.ambig %>% group_by(genomic_loc, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#re-order file 
non_redundant_TEmetadata_low.ambig_sum_genomic_loc_order <- non_redundant_TEmetadata_low.ambig_sum_genomic_loc %>% ungroup(genomic_loc) %>% mutate(genomic_loc = factor(genomic_loc, levels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic"), labels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic")))

#get mean proportion present by bins_2kb of TEs within 10kb of genes 
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin <- non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig %>% group_by(bins_2kb, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#remove genes with bad bin names 
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin <- subset(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin, bins_2kb != "TE_1.5kb_downstream" & bins_2kb != "TE_1kb_upstream" & bins_2kb != "TE_1kb_downstream")
#reorder variables 
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin_order <- non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin %>% ungroup(bins_2kb) %>% mutate(bins_2kb = factor(bins_2kb, levels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb"), labels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb")))

#bin LTR data by prop present 
#bin by prop_present 
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.10] <- 10
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.10 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.20] <- 20
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.20 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.30] <- 30
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.30 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.40] <- 40
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.40 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.50] <- 50
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.50 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.60] <- 60
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.60 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.70] <- 70
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.70 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.80] <- 80
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.80 & non_redundant_TEs_LTRs_low.ambig$prop_present <= 0.90] <- 90
non_redundant_TEs_LTRs_low.ambig$prop_bin[non_redundant_TEs_LTRs_low.ambig$prop_present > 0.90] <- 100
non_redundant_TEs_LTRs_low.ambig_summarize <- non_redundant_TEs_LTRs_low.ambig %>% group_by(prop_bin) %>% summarize(mean_sim = mean(LTR_similarity), std_err = std.error(LTR_similarity), std.dev = sd(LTR_similarity))
non_redundant_TEs_LTRs_low.ambig_summarize_nest <- non_redundant_TEs_LTRs_low.ambig %>% group_by(prop_bin, nested_status) %>% summarize(mean_sim = mean(LTR_similarity), std_err = std.error(LTR_similarity), std.dev = sd(LTR_similarity))

##########################################################################################
##########################################################################################
#GRAPHS: non-redundant dataset
#TE frequency distribution
ggplot(non_redundant_TEmetadata_allTEs_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 508 genotypes by order \nnon-redundant TE dataset") + xlab("Proportion of TEs called as Present") + facet_wrap(. ~ order, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_TEdist_order_scaled.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_TEdist_order_scaled.svg", device = "svg", width = 10, height = 4.5, units = "in")
#prop present nested vs non_nested TEs 
ggplot(non_redundant_TEmetadata_allTEs_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 508 genotypes by order and nested status\nnon-redundant TE dataset") + xlab("Proportion of TEs called as Present") + facet_wrap(nested_status ~ order, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEdist_order_nest_scaled.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEdist_order_nest_scaled.svg", device = "svg", width = 10, height = 6, units = "in")
ggplot(non_redundant_TEmetadata_allTEs_low.ambig, aes(x = prop_present, y = fam_size)) + geom_point(alpha = 0.01) + ggtitle("Non-redundant TEs proportion present vs family size") + facet_grid(order ~ .)
#prop present by genomic location
ggplot(non_redundant_TEmetadata_low.ambig_sum_genomic_loc_order, aes(x = genomic_loc, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by genomic location \n non-redundant TE dataset") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_genomic_loc_wcounts.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_genomic_loc_wcounts.svg", device = "svg", width = 10, height = 6, units = "in")
ggplot(non_redundant_TEmetadata_low.ambig_sum_genomic_loc_order, aes(x = genomic_loc, y = mean_pres)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by genomic location \n non-redundant TE dataset") + facet_grid(order ~ .) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_genomic_loc_nocounts.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_genomic_loc_nocounts.svg", device = "svg", width = 10, height = 6, units = "in")

#graph prop present by genomic location, 2kb bins, only within 10kb of genes 
ggplot(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin_order, aes(x = bins_2kb, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proprotion of genotypes a TE is present in by genomic location \nonly TEs within 10kb of a gene non-redundant TE dataset") + facet_grid(order ~ .) + geom_text(vjust = 2, hjust = 0.2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_gene10kb_wcount.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_gene10kb_wcount.svg", device = "svg", width = 10, height = 7, units = "in")
ggplot(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig_loc_2kb.bin_order, aes(x = bins_2kb, y = mean_pres)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proprotion of genotypes a TE is present in by genomic location \nonly TEs within 10kb of a gene non-redundant TE dataset") + facet_grid(order ~ .) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_gene10kb_nocount.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_gene10kb_nocount.svg", device = "svg", width = 10, height = 7, units = "in")

#LTRs 
ggplot(non_redundant_TEs_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std.dev, ymax = mean_sim + std.dev)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard deviation non-redundant TE dataset") + ylim(87.5,100)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_sd.png", device = "png", width = 7, height = 3.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_sd.svg", device = "svg", width = 7, height = 3.5, units = "in")
ggplot(non_redundant_TEs_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error non-redundant TE dataset") + ylim(87.5,100)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_se.png", device = "png", width = 7, height = 3.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_se.svg", device = "svg", width = 7, height = 3.5, units = "in")
ggplot(non_redundant_TEs_LTRs_low.ambig_summarize_nest, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error non-redundant TE dataset") + ylim(87.5,100) + facet_grid(nested_status ~ .)
ggsave(filename = "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_nest_se.png", device = "png", width = 7, height = 4, units = "in")
ggsave(filename = "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_nest_se.svg", device = "svg", width = 7, height = 4, units = "in")
ggplot(non_redundant_TEs_LTRs_low.ambig_summarize_nest, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std.dev, ymax = mean_sim + std.dev)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error non-redundant TE dataset") + ylim(87.5,100) + facet_grid(nested_status ~ .)
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_nest_sd.png", device = "png", width = 7, height = 4, units = "in")
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_mean_pres_LTR_propbin_nest_sd.svg", device = "svg", width = 7, height = 4, units = "in")

#####################################################
#####################################################
#Correlation between B73 and Mo17 TE numbers 
cor(as.numeric(non_redundant_TEs_rawcounts_B73.Mo17$B73.n_present/(non_redundant_TEs_rawcounts_B73.Mo17$B73.n_present + non_redundant_TEs_rawcounts_B73.Mo17$B73.n_absent)), as.numeric(non_redundant_TEs_rawcounts_B73.Mo17$Mo17.n_present/(non_redundant_TEs_rawcounts_B73.Mo17$Mo17.n_present + non_redundant_TEs_rawcounts_B73.Mo17$Mo17.n_absent)))
#cor = 0.9633707
0.9633707*0.9633707
#Modeling

#Non-redundant TIRs 
nonredundant_TIRs <- subset(non_redundant_TEmetadata_allTEs_low.ambig, order == "TIR")
nonredundant_TIRs_glm <- glm(prop_present ~ as.factor(genomic_loc) + as.factor(nested_status) + TE_len + fam_size, data = nonredundant_TIRs, family = "quasibinomial")
nonredundant_TIRs_nested <- subset(nonredundant_TIRs, nested_status == "nested_TE")
nonredundant_TIRs_nonnested <- subset(nonredundant_TIRs, nested_status == "non_nested_TE")
ks.test(nonredundant_TIRs_nested$prop_present, nonredundant_TIRs_nonnested$prop_present)
#emmeans
library(emmeans)
nonredundant_TIRs_glm.rg <- ref_grid(nonredundant_TIRs_glm)
nonredundant_TIRs_emmeans_loc <- emmeans(nonredundant_TIRs_glm.rg, "genomic_loc", infer = T, level = 0.95)
pairs(nonredundant_TIRs_emmeans_loc)
emmip(nonredundant_TIRs_glm, nested_status ~ genomic_loc) + theme_bw() + labs(y = "Estimated marginal mean \n TE frequency in population", x = "Location in genome") + ggtitle("Estimated marginal means non-redundant TE datset - TIRs") + geom_point(size = 3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_TIRs_emmeans.png", device = "png", width = 8, height = 4.5, units = "in")
plot(nonredundant_TIRs_glm.rg, by = "genomic_loc")
nonredundant_TIRs_emmeans_nest <- emmeans(nonredundant_TIRs_glm.rg, "nested_status", infer = T, level = 0.95)
pairs(nonredundant_TIRs_emmeans_nest)
nonredundant_TIRs_glm_int <- glm(prop_present ~ as.factor(genomic_loc)*nested_status + TE_len*nested_status + fam_size*nested_status, data = nonredundant_TIRs, family = "quasibinomial")
summary(nonredundant_TIRs_glm_int)

nonredundant_Helitrons <- subset(non_redundant_TEmetadata_allTEs_low.ambig, order == "Helitron")
nonredundant_Helitrons_glm <- glm(prop_present ~ genomic_loc + nested_status + TE_len + fam_size, data = nonredundant_Helitrons, family = "quasibinomial")
nonredundant_Helitrons_nested <- subset(nonredundant_Helitrons, nested_status == "nested_TE")
nonredundant_Helitrons_nonnested <- subset(nonredundant_Helitrons, nested_status == "non_nested_TE")
ks.test(nonredundant_Helitrons_nested$prop_present, nonredundant_Helitrons_nonnested$prop_present)
nonredundant_Helitrons_glm.rg <- ref_grid(nonredundant_Helitrons_glm)
nonredundant_Helitrons_emmeans_loc <- emmeans(nonredundant_Helitrons_glm.rg, "genomic_loc", infer = T, level = 0.95)
pairs(nonredundant_Helitrons_emmeans_loc)
nonredundant_Helitrons_emmeans_nest <- emmeans(nonredundant_Helitrons_glm.rg, "nested_status", infer = T, level = 0.95)
pairs(nonredundant_Helitrons_emmeans_nest)
emmip(nonredundant_Helitrons_glm, nested_status ~ genomic_loc) + theme_bw() + labs(y = "Estimated marginal mean \n TE frequency in population", x = "Location in genome") + ggtitle("Estimated marginal means non-redundant TE datset - Helitrons") + geom_point(size = 3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_Helitrons_emmeans.png", device = "png", width = 8, height = 4.5, units = "in")

nonredundant_LTRs_glm <- glm(prop_present ~ genomic_loc + nested_status + TE_len + fam_size + LTR_similarity, data = non_redundant_TEs_LTRs_low.ambig, family = "quasibinomial")
nonredundant_LTRs_nested <- subset(non_redundant_TEs_LTRs_low.ambig, nested_status == "nested_TE")
nonredundant_LTRs_nonnested <- subset(non_redundant_TEs_LTRs_low.ambig, nested_status == "non_nested_TE")
ks.test(nonredundant_LTRs_nested$prop_present, nonredundant_LTRs_nonnested$prop_present)
#GLM: LTRs age* nested 
nonredundant_LTRs_glm_age <- glm(prop_present ~ nested_status * LTR_similarity, data = non_redundant_TEs_LTRs_low.ambig, family = "quasibinomial")
nonredundant_LTRs_glm.rg <- ref_grid(nonredundant_LTRs_glm)
nonredundant_LTRs_emmeans_loc <- emmeans(nonredundant_LTRs_glm.rg, "genomic_loc", infer = T, level = 0.95)
pairs(nonredundant_LTRs_emmeans_loc)
nonredundant_LTRs_emmeans_nest <- emmeans(nonredundant_LTRs_glm.rg, "nested_status", infer = T, level = 0.95)
pairs(nonredundant_LTRs_emmeans_nest)
emmip(nonredundant_LTRs_glm, nested_status ~ genomic_loc) + theme_bw() + labs(y = "Estimated marginal mean \n TE frequency in population", x = "Location in genome") + ggtitle("Estimated marginal means non-redundant TE datset - LTRs") + geom_point(size = 3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename =  "~/Documents/TE_project/non_redundant_TE_set/WiDiv508_nonredundant_set/Figures/non_redundant_TEs_LTRs_emmeans.png", device = "png", width = 8, height = 4.5, units = "in")

#only TEs within 10 kb of genes #there are ~10/20 oddly labeled TEs, will remove "w.in_10kb.downstream" becuase distance to gene is 0 but genomic_loc is TE_5.10 or TE_1.5 kb, only 3 TEs
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig  <- subset(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig, bins_2kb != "w.in_10kb.downstream")
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb[non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb == "TE_1kb_upstream"] <- "up_0kb.2kb"
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb[non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb == "TE_1.5kb_downstream"] <- "down_0kb.2kb"
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb[non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb == "TE_1kb_downstream"] <- "down_0kb.2kb"
non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb[non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig$bins_2kb == "w.in_10kb.upstream"] <- "up_0kb.2kb"

nonredundant_TIRs_10kb <- subset(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig, order == "TIR")
nonredundant_TIRs_10kb_glm <- glm(prop_present ~ as.numeric(distance_to_gene) + bins_2kb + nested_status + TE_len + fam_size, data = nonredundant_TIRs_10kb, family = "quasibinomial")
nonredundant_TIRs_10kb_glm_upneg <- glm(prop_present ~ as.numeric(location_up_neg) + nested_status + TE_len + fam_size, data = nonredundant_TIRs_10kb, family = "quasibinomial")
nonredundant_TIRs_10kb_glm.rg <- ref_grid(nonredundant_TIRs_10kb_glm)
nonredundant_TIRs_10kb_emmeans_loc <- emmeans(nonredundant_TIRs_10kb_glm.rg, "bins_2kb", infer = T, level = 0.95)
pairs(nonredundant_TIRs_10kb_emmeans_loc)

nonredundant_Helitrons_10kb <- subset(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig, order == "Helitron")
nonredundant_Helitrons_10kb_glm <- glm(prop_present ~ as.numeric(distance_to_gene) + bins_2kb + nested_status + TE_len + fam_size, data = nonredundant_Helitrons_10kb, family = "quasibinomial")
nonredundant_Helitrons_10kb_glm_upneg <- glm(prop_present ~ as.numeric(location_up_neg) + nested_status + TE_len + fam_size, data = nonredundant_Helitrons_10kb, family = "quasibinomial")
nonredundant_Helitrons_10kb_glm.rg <- ref_grid(nonredundant_Helitrons_10kb_glm)
nonredundant_Helitrons_10kb_emmeans_loc <- emmeans(nonredundant_Helitrons_10kb_glm.rg, "bins_2kb", infer = T, level = 0.95)
pairs(nonredundant_Helitrons_10kb_emmeans_loc)

nonredundant_LTRs_10kb <- subset(non_redundant_TEmetadata_allTEs_10kb_genes_low.ambig, order == "LTR")
LTR_sim <- data.frame(non_redundant_TEs_LTRs_low.ambig$TE_name, non_redundant_TEs_LTRs_low.ambig$LTR_similarity)
nonredundant_LTRs_10kb <- merge(nonredundant_LTRs_10kb, LTR_sim, by.x = "TE_name", by.y = "non_redundant_TEs_LTRs_low.ambig.TE_name")
nonredundant_LTRs_10kb_glm <- glm(prop_present ~ as.numeric(distance_to_gene) + nested_status + TE_len + bins_2kb + fam_size + non_redundant_TEs_LTRs_low.ambig.LTR_similarity, data = nonredundant_LTRs_10kb, family = "quasibinomial")
nonredundant_LTRs_10kb_glm_upneg <- glm(prop_present ~ as.numeric(location_up_neg) + nested_status + TE_len + fam_size + non_redundant_TEs_LTRs_low.ambig.LTR_similarity, data = nonredundant_LTRs_10kb, family = "quasibinomial")
nonredundant_LTRs_10kb_glm.rg <- ref_grid(nonredundant_LTRs_10kb_glm)
nonredundant_LTRs_10kb_emmeans_loc <- emmeans(nonredundant_LTRs_10kb_glm.rg, "bins_2kb", infer = T, level = 0.95)
pairs(nonredundant_LTRs_10kb_emmeans_loc)


#end of script

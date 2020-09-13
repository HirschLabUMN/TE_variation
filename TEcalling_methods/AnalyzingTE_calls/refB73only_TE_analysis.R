#analysis of TE calling data from 342 WiDiv genotypes - B73 reference 

library(tidyverse)
library(reshape2)
library(ggplot2)
library(scales)
#library(RColorBrewer)
library(dplyr)
library(plotrix)
library(data.table)
#read in data files 
#all these files have high ambig and low ambig TEs 
WiDiv342_refB73_TEmetadata <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_ref.B73_3cat_TE.metadata.txt")
WiDiv342_refB73_TEmetadata_10kb_genes <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_ref.B73_3cat_TE.metadata.bp.togenes_upneg.bins.txt")
WiDiv342_refB73_TEmetadata_LTRs <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_ref.B73_3cat_TE.metadata.LTRsimilarity.txt")
WiDiv342_refB73_TEmetadata_10kb_genes_LTRs <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_ref.B73_3cat_TE.metadata.bp.togenes.LTRsimilarity.txt")

#other data files 
#B73 TEs within 10kb of a gene, with gene name (to connect with gene expression data)
B73_TE_within.10kb.gene.name <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/B73_filteredTE_location_bp.togenes_w.genename_long.ts.nodups.txt")
#gene expression data

#TE expression data - only 10 tissues 
TEfam_exp_select_tissues <- read_tsv("~/Documents/TE_project/gene_te_expression_data/TEfam_select_tissue_exp.txt")
TEfam_exp_select_tissues  <- TEfam_exp_select_tissues %>% rename(exp_count_TEfam = exp_count) %>% rename(exp_prop_TEfam = exp_prop) %>% select(te_fam_name, exp_count_TEfam, exp_prop_TEfam)
B73_gene_exp_select_tissues <- read_tsv("~/Documents/TE_project/gene_te_expression_data/B73_gene_select_tissue_exp.txt")

#combine gene expression data with B73_TE_within.10kb.gene.name file
B73_TE_within.10kb.gene.name_gene.exp <- full_join(B73_TE_within.10kb.gene.name, B73_gene_exp_select_tissues, by = c("gene_name" = "rn"))
#remove rows where TE name is NA and remove rows where distance_to_gene = te_part_ingene
B73_TE_within.10kb.gene.name_gene.exp <- B73_TE_within.10kb.gene.name_gene.exp[!is.na(B73_TE_within.10kb.gene.name_gene.exp$TE_name),]
B73_TE_within.10kb.gene.name_gene.exp <- subset(B73_TE_within.10kb.gene.name_gene.exp, distance_to_gene != "te_part_ingene")
#remove duplicate entries
#use 'distinct'; if there are duplicate rows only the first row is preserved 
#first sort columns so smaller distance_to_gene values show up first
#order sorts assending 
B73_TE_within.10kb.gene.name_gene.exp_uniq <- B73_TE_within.10kb.gene.name_gene.exp[with(B73_TE_within.10kb.gene.name_gene.exp, order(TE_name, distance_to_gene)),] %>% distinct(TE_name, .keep_all = TRUE)
#rename exp_prop and exp_count columns to be distinct from TE information; choose to save only TE_name, gene_name, exp_count and exp_prop columns 
B73_TE_within.10kb.gene.name_gene.exp_uniq <- B73_TE_within.10kb.gene.name_gene.exp_uniq %>% rename(exp_count_gene = exp_count) %>% rename(exp_prop_gene = exp_prop) %>% select(TE_name, gene_name, exp_count_gene, exp_prop_gene)

#add TE expression data to WiDiv metadata files 
WiDiv342_refB73_TEmetadata <- full_join(WiDiv342_refB73_TEmetadata, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv342_refB73_TEmetadata <- WiDiv342_refB73_TEmetadata[!is.na(WiDiv342_refB73_TEmetadata$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv342_refB73_TEmetadata[is.na(WiDiv342_refB73_TEmetadata)] <- 0
#add TE expression data to 3 other data files 
#########################################
WiDiv342_refB73_TEmetadata_10kb_genes <- full_join(WiDiv342_refB73_TEmetadata_10kb_genes, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv342_refB73_TEmetadata_10kb_genes <- WiDiv342_refB73_TEmetadata_10kb_genes[!is.na(WiDiv342_refB73_TEmetadata_10kb_genes$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv342_refB73_TEmetadata_10kb_genes[is.na(WiDiv342_refB73_TEmetadata_10kb_genes)] <- 0
#join mean gene expression values of gene(s) that the TE is closest to 
#use left join: 'return all rows from x, and all columsn from x and y. Rows in x with no match in y will have NA values in the new columns'
#WiDiv342_refB73_TEmetadata_10kb_genes is x and B73_TE_within.10kb.gene.name_gene.exp_sum is y; will replace NA with 0 = any genes close to that TE have no expression 
WiDiv342_refB73_TEmetadata_10kb_genes <- left_join(WiDiv342_refB73_TEmetadata_10kb_genes, B73_TE_within.10kb.gene.name_gene.exp_uniq, by = c("TE_name" = "TE_name"))
WiDiv342_refB73_TEmetadata_10kb_genes[is.na(WiDiv342_refB73_TEmetadata_10kb_genes)] <- 0

#all LTRs
WiDiv342_refB73_TEmetadata_LTRs <- full_join(WiDiv342_refB73_TEmetadata_LTRs, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv342_refB73_TEmetadata_LTRs <- WiDiv342_refB73_TEmetadata_LTRs[!is.na(WiDiv342_refB73_TEmetadata_LTRs$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv342_refB73_TEmetadata_LTRs[is.na(WiDiv342_refB73_TEmetadata_LTRs)] <- 0
#LTRs w/in 10kb of a gene
WiDiv342_refB73_TEmetadata_10kb_genes_LTRs <- full_join(WiDiv342_refB73_TEmetadata_10kb_genes_LTRs, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv342_refB73_TEmetadata_10kb_genes_LTRs <- WiDiv342_refB73_TEmetadata_10kb_genes_LTRs[!is.na(WiDiv342_refB73_TEmetadata_10kb_genes_LTRs$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv342_refB73_TEmetadata_10kb_genes_LTRs[is.na(WiDiv342_refB73_TEmetadata_10kb_genes_LTRs)] <- 0

#choose only low ambig TEs and for TIRs, helitrons and LTRs 
WiDiv342_refB73_TEmetadata_low.ambig <- subset(WiDiv342_refB73_TEmetadata, ambig_cat == "low_ambig" & order == "Helitron" | ambig_cat == "low_ambig" & order == "TIR" | ambig_cat == "low_ambig" & order == "LTR")
WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig <- subset(WiDiv342_refB73_TEmetadata_10kb_genes, ambig_cat == "low_ambig" & order == "Helitron" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "TIR" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "LTR" & distance_to_gene != "te_part_ingene")


#get mean proportion present by genomic location bins of all TEs
WiDiv342_refB73_TEmetadata_low.ambig_sum_genomic_loc <- WiDiv342_refB73_TEmetadata_low.ambig %>% group_by(genomic_loc, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#re-order file 
WiDiv342_refB73_TEmetadata_low.ambig_sum_genomic_loc_order <- WiDiv342_refB73_TEmetadata_low.ambig_sum_genomic_loc %>% ungroup(genomic_loc) %>% mutate(genomic_loc = factor(genomic_loc, levels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic"), labels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic")))
#get mean proportion present by family and then by proportion of tissues a TE family is expressed in 
WiDiv342_refB73_TEmetadata_low.ambig_prop_expressed <- WiDiv342_refB73_TEmetadata_low.ambig %>% group_by(fam_size, family, order, exp_prop_TEfam, exp_count_TEfam) %>% summarize(mean_pres_by_fam = mean(prop_present))

#get mean proportion present by bins_2kb of TEs within 10kb of genes 
WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin <- WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(bins_2kb, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#remove genes with bad bin names 
WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin <- subset(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin, bins_2kb != "TE_1.5kb_downstream" & bins_2kb != "TE_1kb_upstream" & bins_2kb != "TE_1kb_downstream")
#reorder variables 
WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin_order <- WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin %>% ungroup(bins_2kb) %>% mutate(bins_2kb = factor(bins_2kb, levels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb"), labels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb")))
#for TEs within 10kb of genes, look at mean prop present by the percentage of tissues a gene is expressed in 
WiDiv342_refB73_TEmetadata_low.ambig_prop.gene_expressed <- WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_gene, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed <- WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_TEfam, fam_size, family, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())


WiDiv342_refB73_TEmetadata_LTRs_low.ambig <- subset(WiDiv342_refB73_TEmetadata_LTRs, ambig_cat == "low_ambig")
#bin by prop_present 
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.10] <- 10
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.10 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.20] <- 20
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.20 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.30] <- 30
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.30 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.40] <- 40
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.40 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.50] <- 50
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.50 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.60] <- 60
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.60 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.70] <- 70
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.70 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.80] <- 80
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.80 & WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.90] <- 90
WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv342_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.90] <- 100
WiDiv342_refB73_TEmetadata_LTRs_low.ambig_summarize <- WiDiv342_refB73_TEmetadata_LTRs_low.ambig %>% group_by(prop_bin) %>% summarize(mean_sim = mean(LTR_similarity), std_err = std.error(LTR_similarity))

WiDiv342_refB73_TEmetadata_10kb_genes_LTRs_low.ambig <- subset(WiDiv342_refB73_TEmetadata_10kb_genes_LTRs, ambig_cat == "low_ambig")

#graphs 
#TE frequency distribution
ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = prop_present, fill = subgenome)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 342 genotypes \nreference genome = B73") + xlab("Proportion of TEs called as Present") + scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_fillsubgenome.png", device = "png", width = 8, height = 5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_fillsubgenome.svg", device = "svg", width = 8, height = 5, units = "in")
ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 342 genotypes by subgenome \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_wrap(. ~ subgenome, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_subgenome_scaled.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_subgenome_scaled.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 342 genotypes by subgenome \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_grid(. ~ subgenome)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_subgenome.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_subgenome.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 342 genotypes by order \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_grid(. ~ order)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_order.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_order.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 342 genotypes by order \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_wrap(. ~ order, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_order_scaled.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_TEdist_order_scaled.svg", device = "svg", width = 10, height = 4.5, units = "in")
#prop present and prop tissues TE is expressed in
library(viridis)
library(wesanderson)
pal <- wes_palette("Cavalcanti1", n = 10, type = "continuous")
ggplot(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = genomic_loc, y = prop_present, color = exp_count)) + geom_jitter(alpha = 0.02, width = 0.2, height = 0) + ggtitle("TE frequency and distance to gene") + xlab("Distance to gene")  + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_gradient(low = "#b2182b", high = "#2166ac")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_gene10kb_col_expcount.png", device = "png", width = 8, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_gene10kb_col_expcount.svg", device = "svg", width = 8, height = 4.5, units = "in")

#ggplot(WiDiv342_refB73_TEmetadata_low.ambig, aes(x = exp_prop_TEfam, y = prop_present)) + geom_jitter(alpha = 0.02, width = 0.05) + ggtitle("TE frequency and proportion of tissues a TE is expressed in \nall TEs") + xlab("Proportion of tissues a TE family is expressed")  + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Proportion of genotypes a TE is present in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_present_expprop_allTEs.png", device = "png", width = 8, height = 4.5, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_present_expprop_allTEs.svg", device = "svg", width = 8, height = 4.5, units = "in")
#above is incorrect because expression data is by TE family, not individual TE



#graph prop present by genomic location
ggplot(WiDiv342_refB73_TEmetadata_low.ambig_sum_genomic_loc_order, aes(x = genomic_loc, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proprotion of genotypes a TE is present in by genomic location ") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_genomic_loc.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_genomic_loc.svg", device = "svg", width = 10, height = 6, units = "in")
#graph mean proportion present genotypes in each TE family by proportion of tissues a TE family is expressed in
ggplot(WiDiv342_refB73_TEmetadata_low.ambig_prop_expressed, aes(x = exp_prop_TEfam, y = mean_pres_by_fam, color = fam_size)) + geom_jitter(width = 0.05, height = 0) + scale_color_gradient2(low = "#543005", mid = "#f5f5f5", high = "#018571", midpoint = 4000) + ggtitle("Proportion of tissues expressed in a TE family by mean proportion present by TE family \nall TEs") + xlab("Proportion of tissues TE family is expressed in") + ylab("Mean proportion of genotypes a TE is present in by TE family") + facet_grid(order ~ .)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_presentfam_by_propTEfam_exp.png", device = "png", width = 10, height = 8, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_presentfam_by_propTEfam_exp.svg", device = "svg", width = 10, height = 8, units = "in")

#graph prop present by genomic location, 2kb bins, only within 10kb of genes 
ggplot(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin_order, aes(x = bins_2kb, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proprotion of genotypes a TE is present in by genomic location \nonly TEs within 10kb of a gene") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_gene10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_gene10kb.svg", device = "svg", width = 10, height = 7, units = "in")
#graph prop present by prop tissues a TE fam is expressed in, TEs only within 10 kb of genes
ggplot(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed, aes(x = exp_prop_TEfam, y = mean_pres, color = fam_size)) + geom_jitter(alpha = 0.5, width = 0.05, height = 0)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues a TE family is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ .)  + ylim(0,1) + xlab("Proportion of tissues a TE family is expressed in") + scale_color_gradient2(low = "#8c510a", mid = "#f5f5f5", high = "#018571", midpoint = 6000)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_prop_exp_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_prop_exp_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")

ggplot(WiDiv342_refB73_TEmetadata_low.ambig_prop.gene_expressed, aes(x = exp_prop_gene, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.05) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1) + xlab("Proportion of tissues the nearest gene is expressed in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_prop_exp_neargene_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_present_prop_exp_neargene_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")
#graph raw data of propotion of genotypes a TE is in by gene expression percentage
ggplot(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene, y = prop_present)) + geom_jitter(alpha = 0.05, width = 0.02, height = 0) + facet_grid(order ~ .) + ggtitle("Proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nOnlys TEs within 10bk of a gene") 
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_proppresent_exp_neargene_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_proppresent_exp_neargene_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")

#LTRs 
ggplot(WiDiv342_refB73_TEmetadata_LTRs_low.ambig, aes(x = prop_bin, y = LTR_similarity)) + geom_jitter(alpha = 0.02, width = 0.3, height = 0.0) + ggtitle("LTR similarity binned by proportion of genotypes a TE is present in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_prop_present_LTR_propbin.png", device = "png", width = 8, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_prop_present_LTR_propbin.svg", device = "svg", width = 8, height = 4.5, units = "in")
ggplot(WiDiv342_refB73_TEmetadata_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in") + ylim(90,100)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_se_present_LTR_propbin.png", device = "png", width = 8, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv342_analysis/WiDiv342_refB73_mean_se_present_LTR_propbin.svg", device = "svg", width = 8, height = 4.5, units = "in")


###################################################################################################################
###################################################################################################################
#models 
###################################################################################################################
#correlation between distance to gene and population frequency? 
cor(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig$prop_present, as.numeric(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig$distance_to_gene))
#-0.14
cor(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig$prop_present, as.numeric(WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig$fam_size))

#####################################################################################################################
#linear modeling with B73 TEs
#run glm for mix of categorical and numerical predictor variables
#just B73 TIR
B73_TIRs <- subset(WiDiv342_refB73_TEmetadata_low.ambig, order == "TIR")
#add column for 'in gene'
B73_TIRs$genomic_loc_gene <- B73_TIRs$genomic_loc
B73_TIRs$genomic_loc_gene[B73_TIRs$genomic_loc == "5prime_UTR" | B73_TIRs$genomic_loc == "3prime_UTR" | B73_TIRs$genomic_loc == "exon" | B73_TIRs$genomic_loc == "intron" | B73_TIRs$genomic_loc == "TE_encompassed_by_gene" ] <- "in_gene"

B73_TIR_glm <- glm(prop_present ~ factor(genomic_loc_gene) + factor(subgenome) + TE_len + fam_size, data = B73_TIRs)
#just B73 Helitrons
B73_Helitrons <- subset(WiDiv342_refB73_TEmetadata_low.ambig, order == "Helitron")
B73_Helitrons_glm <- glm(prop_present ~ factor(genomic_loc) + factor(subgenome) + TE_len + fam_size, data = B73_TIRs)
#LTRS
B73_LTRs_glm <- glm(prop_present ~ factor(genomic_loc) + factor(subgenome) + TE_len + LTR_similarity + fam_size, data = WiDiv342_refB73_TEmetadata_LTRs_low.ambig)

#model mean prop present by TE family with mean # tissues that family is expressed in
cor(WiDiv342_refB73_TEmetadata_low.ambig_prop_expressed$exp_prop_TEfam, WiDiv342_refB73_TEmetadata_low.ambig_prop_expressed$mean_pres_by_fam)
#not normally distributed, esentially no correlation

#glm for B73 tes witin 10 kb of a gene
#genomic_loc and location_up_neg have the same info in different forms 
#both models below cause R to crash
B73_within10kb_glm <- glm(prop_present ~ factor(subgenome) + TE_len + factor(order) + location_up_neg + fam_size + exp_prop, data = WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig)
#with location_up_neg have upstream as neg numbers and downstream as positive
summary(B73_within10kb_glm)

B73_within10kb_distance.to.gene_glm <- glm(prop_present ~ factor(bins_2kb) * fam_size * exp_prop_gene, data = WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig)
summary(B73_within10kb_distance.to.gene_glm)

#B73 TEs within 10 kb model that doens't crash 
B73_within10kb_distance.to.gene_glm <- glm(prop_present ~ distance_to_gene, data = WiDiv342_refB73_TEmetadata_10kb_genes_low.ambig)
#end of script

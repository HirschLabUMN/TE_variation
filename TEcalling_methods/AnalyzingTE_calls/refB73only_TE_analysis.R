#analysis of TE calling data from 508 WiDiv genotypes - B73 reference 

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
WiDiv508_refB73_TEmetadata <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_ref.B73_3cat_TE.metadata.txt")
WiDiv508_refB73_TEmetadata_10kb_genes <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_ref.B73_3cat_TE.metadata.bp.togenes_upneg.bins.txt")
WiDiv508_refB73_TEmetadata_LTRs <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_ref.B73_3cat_TE.metadata.LTRsimilarity.txt")
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_ref.B73_3cat_TE.metadata.bp.togenes.LTRsimilarity.txt")

#other data files 
#B73 TEs witWiDiv508hin 10kb of a gene, with gene name (to connect with gene expression data)
B73_TE_within.10kb.gene.name <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/ref_B73_TEcalling/WiDiv508_analysis/B73_filteredTE_location_bp.togenes_w.genename_long.ts.nodups_psuedoUTRs.txt")
#gene expression data

#TE expression data - only 10 tissues 
TEfam_exp_select_tissues <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/gene_te_expression_data/TEfam_select_tissue_exp.txt")
TEfam_exp_select_tissues  <- TEfam_exp_select_tissues %>% rename(exp_count_TEfam = exp_count) %>% rename(exp_prop_TEfam = exp_prop) %>% select(te_fam_name, exp_count_TEfam, exp_prop_TEfam)
B73_gene_exp_select_tissues <- read_tsv("~/Dropbox/HirschLab_MaizeWork/TE_project/gene_te_expression_data/B73_gene_select_tissue_exp.txt")

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
WiDiv508_refB73_TEmetadata <- full_join(WiDiv508_refB73_TEmetadata, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv508_refB73_TEmetadata <- WiDiv508_refB73_TEmetadata[!is.na(WiDiv508_refB73_TEmetadata$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv508_refB73_TEmetadata[is.na(WiDiv508_refB73_TEmetadata)] <- 0
#add TE expression data to 3 other data files 
#########################################
WiDiv508_refB73_TEmetadata_10kb_genes <- full_join(WiDiv508_refB73_TEmetadata_10kb_genes, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv508_refB73_TEmetadata_10kb_genes <- WiDiv508_refB73_TEmetadata_10kb_genes[!is.na(WiDiv508_refB73_TEmetadata_10kb_genes$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv508_refB73_TEmetadata_10kb_genes[is.na(WiDiv508_refB73_TEmetadata_10kb_genes)] <- 0
#join mean gene expression values of gene(s) that the TE is closest to 
#use left join: 'return all rows from x, and all columsn from x and y. Rows in x with no match in y will have NA values in the new columns'
#WiDiv508_refB73_TEmetadata_10kb_genes is x and B73_TE_within.10kb.gene.name_gene.exp_sum is y; will replace NA with 0 = any genes close to that TE have no expression 
WiDiv508_refB73_TEmetadata_10kb_genes <- left_join(WiDiv508_refB73_TEmetadata_10kb_genes, B73_TE_within.10kb.gene.name_gene.exp_uniq, by = c("TE_name" = "TE_name"))
WiDiv508_refB73_TEmetadata_10kb_genes[is.na(WiDiv508_refB73_TEmetadata_10kb_genes)] <- 0

#all LTRs
WiDiv508_refB73_TEmetadata_LTRs <- full_join(WiDiv508_refB73_TEmetadata_LTRs, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv508_refB73_TEmetadata_LTRs <- WiDiv508_refB73_TEmetadata_LTRs[!is.na(WiDiv508_refB73_TEmetadata_LTRs$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv508_refB73_TEmetadata_LTRs[is.na(WiDiv508_refB73_TEmetadata_LTRs)] <- 0
#LTRs w/in 10kb of a gene
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs <- full_join(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs, TEfam_exp_select_tissues, by = c("family" = "te_fam_name"))
#remove rows where TE name is NA
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs <- WiDiv508_refB73_TEmetadata_10kb_genes_LTRs[!is.na(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$TE_name),]
#replace all 'NA' with 0 - NA are TE families that don't have any expression data 
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs[is.na(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs)] <- 0

#choose only low ambig TEs and for TIRs, helitrons and LTRs 
WiDiv508_refB73_TEmetadata_low.ambig <- subset(WiDiv508_refB73_TEmetadata, ambig_cat == "low_ambig" & order == "Helitron" | ambig_cat == "low_ambig" & order == "TIR" | ambig_cat == "low_ambig" & order == "LTR")
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig <- subset(WiDiv508_refB73_TEmetadata_10kb_genes, ambig_cat == "low_ambig" & order == "Helitron" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "TIR" & distance_to_gene != "te_part_ingene"| ambig_cat == "low_ambig" & order == "LTR" & distance_to_gene != "te_part_ingene")
#methylation data - choose side with higher methylation
WiDiv508_refB73_TEmetadata_low.ambig$cg_high <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_up > WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_up, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down)
WiDiv508_refB73_TEmetadata_low.ambig$chh_high <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_up > WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_up, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down)
WiDiv508_refB73_TEmetadata_low.ambig$chg_high <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_up > WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_down, WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_up, WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_down)
#the same, but with TEs within 10 kb of genes 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$cg_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_up > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_down, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_up, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_down)
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$chh_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_up > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_down, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_up, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_down)
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$chg_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_up > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_down, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_up, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_down)
#methylation data for LTRs 
WiDiv508_refB73_TEmetadata_LTRs$cg_high <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_cg_up > WiDiv508_refB73_TEmetadata_LTRs$mean_cg_down, WiDiv508_refB73_TEmetadata_LTRs$mean_cg_up, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down)
WiDiv508_refB73_TEmetadata_LTRs$chh_high <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_chh_up > WiDiv508_refB73_TEmetadata_LTRs$mean_chh_down, WiDiv508_refB73_TEmetadata_LTRs$mean_chh_up, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down)
WiDiv508_refB73_TEmetadata_LTRs$chg_high <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_chg_up > WiDiv508_refB73_TEmetadata_LTRs$mean_chg_down, WiDiv508_refB73_TEmetadata_LTRs$mean_chg_up, WiDiv508_refB73_TEmetadata_LTRs$mean_chg_down)
#LTRs only within 10kb of genes 
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$cg_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_up > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_down, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_up, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_down)
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$chh_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_up > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_down, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_up, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_down)
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$chg_high <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_up > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_down, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_up, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_down)

#methylation data - choose side with higher methylation, mean methylation by family 
WiDiv508_refB73_TEmetadata_low.ambig$cg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_up_fam > WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_up_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down_fam)
WiDiv508_refB73_TEmetadata_low.ambig$chh_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_up_fam > WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_up_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down_fam)
WiDiv508_refB73_TEmetadata_low.ambig$chg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_up_fam > WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_down_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_up_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_chg_down_fam)
#the same, but with TEs within 10 kb of genes 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$cg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_cg_down_fam)
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$chh_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chh_down_fam)
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$chg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$mean_chg_down_fam)
#methylation data for LTRs 
WiDiv508_refB73_TEmetadata_LTRs$cg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_cg_up_fam > WiDiv508_refB73_TEmetadata_LTRs$mean_cg_down_fam, WiDiv508_refB73_TEmetadata_LTRs$mean_cg_up_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_cg_down_fam)
WiDiv508_refB73_TEmetadata_LTRs$chh_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_chh_up_fam > WiDiv508_refB73_TEmetadata_LTRs$mean_chh_down_fam, WiDiv508_refB73_TEmetadata_LTRs$mean_chh_up_fam, WiDiv508_refB73_TEmetadata_low.ambig$mean_chh_down_fam)
WiDiv508_refB73_TEmetadata_LTRs$chg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_LTRs$mean_chg_up_fam > WiDiv508_refB73_TEmetadata_LTRs$mean_chg_down_fam, WiDiv508_refB73_TEmetadata_LTRs$mean_chg_up_fam, WiDiv508_refB73_TEmetadata_LTRs$mean_chg_down_fam)
#LTRs only within 10kb of genes 
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$cg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_cg_down_fam)
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$chh_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chh_down_fam)
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$chg_high_fam <- ifelse(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_up_fam > WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_down_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_up_fam, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs$mean_chg_down_fam)


#get mean proportion present by genomic location bins of all TEs
WiDiv508_refB73_TEmetadata_low.ambig_sum_genomic_loc <- WiDiv508_refB73_TEmetadata_low.ambig %>% group_by(genomic_loc, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#re-order file 
WiDiv508_refB73_TEmetadata_low.ambig_sum_genomic_loc_order <- WiDiv508_refB73_TEmetadata_low.ambig_sum_genomic_loc %>% ungroup(genomic_loc) %>% mutate(genomic_loc = factor(genomic_loc, levels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic"), labels = c("TE_10.5kb_upstream", "TE_5.1kb_upstream", "TE_1kb_upstream", "5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "TE_1kb_downstream", "TE_1.5kb_downstream", "TE_5.10kb_downstream", "TE_intergenic")))
#get mean proportion present by family and then by proportion of tissues a TE family is expressed in 
WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed <- WiDiv508_refB73_TEmetadata_low.ambig %>% group_by(fam_size, family, order, exp_prop_TEfam, exp_count_TEfam) %>% summarize(mean_pres_by_fam = mean(prop_present))

#get mean proportion present by bins_2kb of TEs within 10kb of genes 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(bins_2kb, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
#remove genes with bad bin names 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin <- subset(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin, bins_2kb != "TE_1.5kb_downstream" & bins_2kb != "TE_1kb_upstream" & bins_2kb != "TE_1kb_downstream")
#reorder variables 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin_order <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin %>% ungroup(bins_2kb) %>% mutate(bins_2kb = factor(bins_2kb, levels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb"), labels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb")))
#for TEs within 10kb of genes, look at mean prop present by the percentage of tissues a gene is expressed in 
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_gene, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_wlocation <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_gene, order, bins_2kb) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_wlocation_order <- WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_wlocation %>% ungroup(bins_2kb) %>% mutate(bins_2kb = factor(bins_2kb, levels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb"), labels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb")))
#create 5 bins for prop tissues a gene is expressed in: 0-2, 3-4, 5-6, 7-8, 9-10 and get means 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene_bin[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene <= 0.2] <- 0.2
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene_bin[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene > 0.2 & WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene <= 0.4 ] <- 0.4
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene_bin[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene > 0.4 & WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene <= 0.6 ] <- 0.6
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene_bin[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene > 0.6 & WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene <= 0.8 ] <- 0.8
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene_bin[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$exp_prop_gene > 0.8] <- 1
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_gene_bin, order, bins_2kb) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order <- WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation %>% ungroup(bins_2kb) %>% mutate(bins_2kb = factor(bins_2kb, levels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb"), labels = c("up_8kb.10kb", "up_6kb.8kb", "up_4kb.6kb", "up_2kb.4kb", "up_0kb.2kb","5prime_UTR", "intron", "exon", "3prime_UTR", "TE_encompassed_by_gene", "TE_encompassing_gene", "down_0kb.2kb", "down_2kb.4kb", "down_4kb.6kb", "down_6kb.8kb", "down_8kb.10kb")))

#mean present by TE family tissue expression
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig %>% group_by(exp_prop_TEfam, fam_size, family, order) %>% summarize(mean_pres = mean(prop_present), std_err = std.error(prop_present), count = n())

WiDiv508_refB73_TEmetadata_LTRs_low.ambig <- subset(WiDiv508_refB73_TEmetadata_LTRs, ambig_cat == "low_ambig")
#bin by prop_present 
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.10] <- 10
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.10 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.20] <- 20
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.20 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.30] <- 30
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.30 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.40] <- 40
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.40 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.50] <- 50
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.50 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.60] <- 60
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.60 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.70] <- 70
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.70 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.80] <- 80
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.80 & WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present <= 0.90] <- 90
WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_bin[WiDiv508_refB73_TEmetadata_LTRs_low.ambig$prop_present > 0.90] <- 100
WiDiv508_refB73_TEmetadata_LTRs_low.ambig_summarize <- WiDiv508_refB73_TEmetadata_LTRs_low.ambig %>% group_by(prop_bin) %>% summarize(mean_sim = mean(LTR_similarity), std_err = std.error(LTR_similarity), std.dev = sd(LTR_similarity))

WiDiv508_refB73_TEmetadata_LTRs_nest_low.ambig_summarize <- WiDiv508_refB73_TEmetadata_LTRs_low.ambig %>% group_by(prop_bin, nested_status) %>% summarize(mean_sim = mean(LTR_similarity), std_err = std.error(LTR_similarity), std.dev = sd(LTR_similarity))

WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig <- subset(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs, ambig_cat == "low_ambig")

##########################################################################################
##########################################################################################
#GRAPHS 
#TE frequency distribution
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present, fill = subgenome)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 509 genotypes \nreference genome = B73") + xlab("Proportion of TEs called as Present") + scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_fillsubgenome.png", device = "png", width = 8, height = 5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_fillsubgenome.svg", device = "svg", width = 8, height = 5, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 509 genotypes by subgenome \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_wrap(. ~ subgenome, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_subgenome_scaled.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_subgenome_scaled.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 509 genotypes by subgenome \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_grid(. ~ subgenome)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_subgenome.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_subgenome.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 509 genotypes by order \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_grid(. ~ order)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order.svg", device = "svg", width = 10, height = 4.5, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 508 genotypes by order \nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_wrap(. ~ order, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order_scaled.png", device = "png", width = 10, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order_scaled.svg", device = "svg", width = 10, height = 4.5, units = "in")
#prop present nested vs non_nested TEs 
ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = prop_present)) + geom_histogram(bins = 50) + ggtitle("TE frequency distribution of TEs called in 508 genotypes by order and nested status\nreference genome = B73") + xlab("Proportion of TEs called as Present") + facet_wrap(nested_status ~ order, scales = "free_y")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order_nest_scaled.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEdist_order_nest_scaled.svg", device = "svg", width = 10, height = 6, units = "in")

#prop present and prop tissues TE is expressed in
library(viridis)
library(wesanderson)
pal <- wes_palette("Cavalcanti1", n = 10, type = "continuous")
ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = genomic_loc, y = prop_present, color = exp_count_gene)) + geom_jitter(alpha = 0.02, width = 0.2, height = 0) + ggtitle("TE frequency and distance to gene") + xlab("Distance to gene")  + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_gradient(low = "#b2182b", high = "#2166ac")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_gene10kb_col_expcount.png", device = "png", width = 8, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_gene10kb_col_expcount.svg", device = "svg", width = 8, height = 4.5, units = "in")

#ggplot(WiDiv508_refB73_TEmetadata_low.ambig, aes(x = exp_prop_TEfam, y = prop_present)) + geom_jitter(alpha = 0.02, width = 0.05) + ggtitle("TE frequency and proportion of tissues a TE is expressed in \nall TEs") + xlab("Proportion of tissues a TE family is expressed")  + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Proportion of genotypes a TE is present in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_present_expprop_allTEs.png", device = "png", width = 8, height = 4.5, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_present_expprop_allTEs.svg", device = "svg", width = 8, height = 4.5, units = "in")
#above is incorrect because expression data is by TE family, not individual TE

#graph prop present by genomic location
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_sum_genomic_loc_order, aes(x = genomic_loc, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by genomic location ") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_genomic_loc.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_genomic_loc.svg", device = "svg", width = 10, height = 6, units = "in")
#graph mean proportion present genotypes in each TE family by proportion of tissues a TE family is expressed in
#add new column for family size, anything over 10k, set to 10k 
WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed$fam_size_limit <- WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed$fam_size
WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed$fam_size_limit[WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed$fam_size > 10000] <- 10000

#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed, aes(x = exp_prop_TEfam, y = mean_pres_by_fam, color = fam_size)) + geom_jitter(width = 0.05, height = 0) + scale_color_gradient2(low = "#543005", mid = "#f5f5f5", high = "#018571", midpoint = 4000) + ggtitle("Proportion of tissues expressed in a TE family by mean proportion present by TE family \nall TEs") + xlab("Proportion of tissues TE family is expressed in") + ylab("Mean proportion of genotypes a TE is present in by TE family") + facet_grid(order ~ .)
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfam_exp.png", device = "png", width = 10, height = 8, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfam_exp.svg", device = "svg", width = 10, height = 8, units = "in")
#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed, aes(x = exp_prop_TEfam, y = mean_pres_by_fam, color = fam_size_limit, size = fam_size_limit)) + geom_jitter(width = 0.05, height = 0, alpha = 0.75) + scale_color_gradient2(low = "#c7eae5", mid = "#35978f", high = "#003c30", midpoint = 5000) + ggtitle("Proportion of tissues expressed in a TE family by mean proportion present by TE family \nall TEs") + xlab("Proportion of tissues TE family is expressed in") + ylab("Mean proportion of genotypes a TE is present in by TE family") + facet_grid(order ~ .)
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_cs.png", device = "png", width = 10, height = 8, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_cs.svg", device = "svg", width = 10, height = 8, units = "in")
#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed, aes(x = exp_prop_TEfam, y = mean_pres_by_fam, color = fam_size_limit)) + geom_jitter(width = 0.05, height = 0) + scale_color_gradient2(low = "#c7eae5", mid = "#35978f", high = "#003c30", midpoint = 5000) + ggtitle("Proportion of tissues expressed in a TE family by mean proportion present by TE family \nall TEs") + xlab("Proportion of tissues TE family is expressed in") + ylab("Mean proportion of genotypes a TE is present in by TE family") + facet_grid(order ~ .)
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_c.png", device = "png", width = 10, height = 8, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_c.svg", device = "svg", width = 10, height = 8, units = "in")
#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop_expressed, aes(x = exp_prop_TEfam, y = mean_pres_by_fam, size = fam_size_limit)) + geom_jitter(width = 0.05, height = 0, alpha = 0.2) + scale_color_gradient2(low = "#c7eae5", mid = "#35978f", high = "#003c30", midpoint = 5000) + ggtitle("Proportion of tissues expressed in a TE family by mean proportion present by TE family \nall TEs") + xlab("Proportion of tissues TE family is expressed in") + ylab("Mean proportion of genotypes a TE is present in by TE family") + facet_grid(order ~ .)
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_s.png", device = "png", width = 10, height = 8, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_presentfam_by_propTEfamexp_s.svg", device = "svg", width = 10, height = 8, units = "in")

#graph prop present by genomic location, 2kb bins, only within 10kb of genes 
ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_sum_loc_2kb.bin_order, aes(x = bins_2kb, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proprotion of genotypes a TE is present in by genomic location \nonly TEs within 10kb of a gene") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1.1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_gene10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_gene10kb.svg", device = "svg", width = 10, height = 7, units = "in")

#graph prop present by prop tissues a TE fam is expressed in, TEs only within 10 kb of genes
#WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed$fam_size_limit <- WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed$fam_size
#WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed$fam_size_limit[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed$fam_size > 10000] <- 10000
#ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig_prop.TEfam_expressed, aes(x = exp_prop_TEfam, y = mean_pres, color = fam_size_limit, size = fam_size_limit)) + geom_jitter(alpha = 0.5, width = 0.05, height = 0)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues a TE family is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ .)  + ylim(0,1) + xlab("Proportion of tissues a TE family is expressed in") + scale_color_gradient2(low = "#c7eae5", mid = "#35978f", high = "#003c30", midpoint = 5000)
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")

ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed, aes(x = exp_prop_gene, y = mean_pres, label = count)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.05) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ .) + geom_text(vjust = 2) + ylim(0,1) + xlab("Proportion of tissues the nearest gene is expressed in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")
#graph the mean proportion of genotyes a TE is expressed in by proportion of tissues the nearest gene is expressed in 
#WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order <- na.omit(WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order)
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order, aes(x = exp_prop_gene_bin, y = mean_pres)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.05) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ bins_2kb)  + ylim(0,1) + xlab("Proportion of tissues the nearest gene is expressed in") + ylab("Mean proportion of genotypes a TE is present in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb_bin2kb.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb_bin2kb.svg", device = "svg", width = 12, height = 7, units = "in")
#same as above, but only 0-2 kb and within genes 
WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order_genes <- subset(WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order, bins_2kb == "up_0kb.2kb" | bins_2kb == "5prime_UTR" | bins_2kb == "intron" | bins_2kb == "exon" | bins_2kb == "3prime_UTR" | bins_2kb == "TE_encompassed_by_gene" | bins_2kb == "TE_encompassing_gene" | bins_2kb == "down_0kb.2kb")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_prop.gene_expressed_tissue.bin_wlocation_order_genes, aes(x = exp_prop_gene_bin, y = mean_pres)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_pres-std_err, ymax = mean_pres + std_err), width = 0.05) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Mean proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nonly TEs within 10 kb of a gene") + facet_grid(order ~ bins_2kb)  + ylim(0,1) + xlab("Proportion of tissues the nearest gene is expressed in") + ylab("Mean proportion of genotypes a TE is present in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb_innext_genes.png", device = "png", width = 10, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_present_prop_exp_neargene_TEs_win10kb_innext_genes.svg", device = "svg", width = 10, height = 6, units = "in")

#graph raw data of propotion of genotypes a TE is in by gene expression percentage
#ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene, y = prop_present)) + geom_jitter(alpha = 0.05, width = 0.02, height = 0) + facet_grid(order ~ .) + ggtitle("Proportion of genotypes a TE is present in by proportion of tissues the nearest gene is expressed in \nOnlys TEs within 10bk of a gene") 
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_exp_neargene_TEs_win10kb.png", device = "png", width = 10, height = 7, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_exp_neargene_TEs_win10kb.svg", device = "svg", width = 10, height = 7, units = "in")

#methylation data 
#only graph non-nested TEs 
WiDiv508_refB73_TEmetadata_low.ambig_non_nest <- subset(WiDiv508_refB73_TEmetadata_low.ambig, nested_status == "non_nested_TE")
#bin by prop_present 
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.10] <- 10
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.10 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.20] <- 20
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.20 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.30] <- 30
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.30 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.40] <- 40
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.40 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.50] <- 50
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.50 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.60] <- 60
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.60 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.70] <- 70
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.70 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.80] <- 80
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.80 & WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present <= 0.90] <- 90
WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_bin[WiDiv508_refB73_TEmetadata_low.ambig_non_nest$prop_present > 0.90] <- 100
WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize <- WiDiv508_refB73_TEmetadata_low.ambig_non_nest %>% group_by(prop_bin, order) %>% summarize(mean_chh = mean(chh_high), std_err_chh = std.error(chh_high), std_dev_chh = sd(chh_high), mean_chg = mean(chg_high), std_err_chg = std.error(chg_high), std_dev_chg = sd(chg_high), mean_cg = mean(cg_high), std_err_cg = std.error(cg_high), std_dev_cg = sd(cg_high), mean_chh_fam = mean(chh_high_fam), std_err_chh_fam = std.error(chh_high_fam), std_dev_chh_fam = sd(chh_high_fam), mean_chg_fam = mean(chg_high_fam), std_err_chg_fam = std.error(chg_high_fam), std_dev_chg_fam = sd(chg_high_fam), mean_cg_fam = mean(cg_high_fam), std_err_cg_fam = std.error(cg_high_fam), std_dev_cg_fam = sd(cg_high_fam))

#ggplot(WiDiv508_refB73_TEmetadata_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error") + ylim(87.5,100)

ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_chh)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_chh - std_dev_chh, ymax = mean_chh + std_dev_chh)) + facet_grid(order ~ .) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\n non-nested TEs only") + ylab("CHH methylation ratio") + xlab("Proportion of genotypes the TE is present in") + ylim(-0.05,0.2)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHH_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHHmethylation.svg", device = "svg", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_chg)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_chg - std_dev_chg, ymax = mean_chg + std_dev_chg)) + facet_grid(order ~ .) + ggtitle("Mean CHG methylation over 1kb up or downstream (whichever side is higher) of the TE\n non-nested TEs only") + ylab("CHG methylation ratio") + xlab("Proportion of genotypes the TE is present in") + ylim(-0.15,1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHG_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHGmethylation.svg", device = "svg", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_cg)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_cg - std_dev_cg, ymax = mean_cg + std_dev_cg)) + facet_grid(order ~ .) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in") 
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CG_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CGmethylation.svg", device = "svg", width = 12, height = 7, units = "in")

ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_cg)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_cg - std_err_cg, ymax = mean_cg + std_err_cg)) + facet_grid(order ~ .) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in") 

#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest, aes(x = prop_present, y = chh_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE \n mean methylation by TE family; non-nested TEs only") + ylab("CHH methylation ratio") + xlab("Proportion of genotypes the TE is present in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHHmethylation_fam_loc.png", device = "png", width = 12, height = 7, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHHmethylation_fam_loc.svg", device = "svg", width = 12, height = 7, units = "in")
#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest, aes(x = prop_present, y = chg_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHG methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CHG methylation ratio") + xlab("Proportion of genotypes the TE is present in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHGmethylation_fam_loc.png", device = "png", width = 12, height = 7, units = "in")
#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHGmethylation_fam_loc.svg", device = "svg", width = 12, height = 7, units = "in")

#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest, aes(x = prop_present, y = cg_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in")
#ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_cg_fam)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_cg_fam - std_err_cg_fam, ymax = mean_cg_fam + std_err_cg_fam)) + facet_grid(order ~ .) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in") 

#ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CGmethylation_fam_loc.png", device = "png", width = 12, height = 7, units = "in")
#ggsave(filename = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CGmethylation_fam_loc.svg", device = "svg", width = 12, height = 7, units = "in")

ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_chh_fam)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_chh_fam - std_dev_chh_fam, ymax = mean_chh_fam + std_dev_chh_fam)) + facet_grid(order ~ .) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in") 
+ ylim(0,1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHH_fam_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHHmethylation_fam.svg", device = "svg", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_chg_fam)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_chg_fam - std_dev_chg_fam, ymax = mean_chg_fam + std_dev_chg_fam)) + facet_grid(order ~ .) + ggtitle("Mean CHG methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CHG methylation ratio") + xlab("Proportion of genotypes the TE is present in") + ylim(0,1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHG_fam_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CHGmethylation_fam.svg", device = "svg", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest, aes(x = prop_present, y = cg_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ .) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in")
ggplot(WiDiv508_refB73_TEmetadata_low.ambig_non_nest_summarize, aes(x = prop_bin, y = mean_cg_fam)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_cg_fam - std_dev_cg_fam, ymax = mean_cg_fam + std_dev_cg_fam)) + facet_grid(order ~ .) + ggtitle("Mean CG methylation over 1kb up or downstream (whichever side is higher) of the TE\n mean methylation by TE family; non-nested TEs only") + ylab("CG methylation ratio") + xlab("Proportion of genotypes the TE is present in") + ylim(0.25,1)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CG_fam_meanprop.png", device = "png", width = 12, height = 7, units = "in")
ggsave(filename = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_proppresent_CG_fam_meanprop.svg", device = "svg", width = 12, height = 7, units = "in")


#methylation and the number of tissues the nearest gene is expressed in
#non nested TEs only
ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene_bin, y = chh_high)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\noly TEs within 10kb of genes") + ylab("CHH methylation ratio") + xlab("Proportion of tissues the nearest gene is expressed in") + ylim(0,0.75)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_propgeneexp_CHHmethylation.png", device = "png", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene_bin, y = chg_high)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\noly TEs within 10kb of genes") + ylab("CHG methylation ratio") + xlab("Proportion of tissues the nearest gene is expressed in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_propgeneexp_CHGmethylation.png", device = "png", width = 12, height = 7, units = "in")

ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene_bin, y = chh_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\noly TEs within 10kb of genes") + ylab("CHH methylation ratio") + xlab("Proportion of tissues the nearest gene is expressed in") + ylim(0,0.5)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_propgeneexp_CHHmethylation.png", device = "png", width = 12, height = 7, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, aes(x = exp_prop_gene_bin, y = chg_high_fam)) + geom_point(alpha = 0.1) + facet_grid(order ~ genomic_loc) + ggtitle("Mean CHH methylation over 1kb up or downstream (whichever side is higher) of the TE\noly TEs within 10kb of genes") + ylab("CHG methylation ratio") + xlab("Proportion of tissues the nearest gene is expressed in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_propgeneexp_CHGmethylation.png", device = "png", width = 12, height = 7, units = "in")

#LTRs 
ggplot(WiDiv508_refB73_TEmetadata_LTRs_low.ambig, aes(x = prop_bin, y = LTR_similarity)) + geom_jitter(alpha = 0.02, width = 0.3, height = 0.0) + ggtitle("LTR similarity binned by proportion of genotypes a TE is present in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_prop_present_LTR_propbin.png", device = "png", width = 8, height = 4.5, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_prop_present_LTR_propbin.svg", device = "svg", width = 8, height = 4.5, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std.dev, ymax = mean_sim + std.dev)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard deviation") + ylim(87.5,100)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_se_present_LTR_propbin_sd.png", device = "png", width = 6, height = 4, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_se_present_LTR_propbin_sd.svg", device = "svg", width = 6, height = 4, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_LTRs_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error") + ylim(87.5,100)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_se_present_LTR_propbin_se.png", device = "png", width = 6, height = 4, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_mean_se_present_LTR_propbin_se.svg", device = "svg", width = 6, height = 4, units = "in")

#nested vs non-nested LTRs 
ggplot(WiDiv508_refB73_TEmetadata_LTRs_nest_low.ambig_summarize, aes(x = prop_bin, y = mean_sim))+ 
  facet_grid(nested_status ~ .)+ geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std.dev, ymax = mean_sim + std.dev)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard deviation") + ylim(87.5,101)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEs_mean_pres_LTR_propbin_nest_sd.png", device = "png", width = 7, height = 4, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEs_mean_pres_LTR_propbin_nest_sd.svg", device = "svg", width = 7, height = 4, units = "in")
ggplot(WiDiv508_refB73_TEmetadata_LTRs_nest_low.ambig_summarize, aes(x = prop_bin, y = mean_sim)) + geom_point(size = 2) + geom_errorbar(aes(ymin = mean_sim - std_err, ymax = mean_sim + std_err)) + ggtitle("Mean LTR similarity binned by proportion of genotypes a TE is present in\nmean +/- standard error") + ylim(87.5,100) + facet_grid(nested_status ~ .)
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEs_mean_pres_LTR_propbin_nest_se.png", device = "png", width = 7, height = 4, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEs_mean_pres_LTR_propbin_nest_se.svg", device = "svg", width = 7, height = 4, units = "in")



#nested vs non nested TEs 
WiDiv509_nested_v_outer_present <- read_tsv("~/Documents/TE_project/ref_B73_TEcalling/WiDiv509_B73_nested.v.outer_present.txt")
WiDiv509_nested_v_outer_present_lowambig <- subset(WiDiv509_nested_v_outer_present, nestedTE_ambigcat == "low_ambig" & outerTE_ambigcat == "low_ambig")
ggplot(WiDiv509_nested_v_outer_present_lowambig, aes(x = nestedTE_prop_present, y = outerTE_prop_present)) + geom_point(alpha = 0.02) + ggtitle("Proportion of genotypes a TE is present in \nnested TEs vs the TE they are nested in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv509_refB73_nested_v_outer.png", device = "png", width = 7, height = 8, units = "in")
#re-do with WiDiv 508 numbers 
WiDiv508_prop_only <- WiDiv508_refB73_TEmetadata[,c(1,20,21)]
WiDiv508_nested_v_outer_int <- merge(WiDiv509_nested_v_outer_present, WiDiv508_prop_only, by.x = "nestedTE_name", by.y = "TE_name")
WiDiv508_nested_v_outer_int <- WiDiv508_nested_v_outer_int[,c(1,4,7,8)]
WiDiv508_nested_v_outer_int <- WiDiv508_nested_v_outer_int %>% rename(nestedTE_prop_present  = prop_present, nestedTE_ambigcat = ambig_cat)
WiDiv508_nested_v_outer <- merge(WiDiv508_nested_v_outer_int, WiDiv508_prop_only, by.x = "outerTE_name", by.y = "TE_name")
WiDiv508_nested_v_outer <- WiDiv508_nested_v_outer %>% rename(outerTE_prop_present = prop_present, outerTE_ambigcat = ambig_cat)
WiDiv508_nested_v_outer_present_lowambig <- subset(WiDiv508_nested_v_outer, nestedTE_ambigcat == "low_ambig" & outerTE_ambigcat == "low_ambig")
ggplot(WiDiv508_nested_v_outer_present_lowambig, aes(x = nestedTE_prop_present, y = outerTE_prop_present)) + geom_point(alpha = 0.02) + ggtitle("Proportion of genotypes a TE is present in B73 reference genome")+
  xlab("Proportion of genotypes the nested TE is Present in") + ylab("Proportion of genotypes the outer TE is Present in") + theme_bw()
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_nested_v_outer.png", device = "png", width = 6, height = 6, units = "in")
ggsave(filename =  "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_nested_v_outer.svg", device = "svg", width = 6, height = 6, units = "in")

cor(WiDiv508_nested_v_outer_present_lowambig$nestedTE_prop_present, WiDiv508_nested_v_outer_present_lowambig$outerTE_prop_present)
0.7599042*0.7599042
####################################################################################################################
####################################################################################################################
#models 
####################################################################################################################
####################################################################################################################
#linear modeling with B73 TEs
#run glm for mix of categorical and numerical predictor variables
#should use family = "binomial" 
#Why? https://stats.stackexchange.com/questions/284843/percentage-as-dependent-variable-in-multiple-linear-regression
#binomial give errors 
#https://stackoverflow.com/questions/12953045/warning-non-integer-successes-in-a-binomial-glm-survey-packages
#can use quasibinomial instead
#using quasibinomial doesn't change results but does get rid of errors 

#just B73 TIR
B73_TIRs <- subset(WiDiv508_refB73_TEmetadata_low.ambig, order == "TIR")
#add column for 'in gene'
B73_TIR_glm <- glm(prop_present ~ genomic_loc + subgenome + nested_status + TE_len + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_TIRs, family = "quasibinomial")
B73_TIR_glm_methlation <- glm(prop_present ~ genomic_loc + subgenome + nested_status + TE_len + fam_size + cg_high + chh_high + chg_high, data = B73_TIRs, family = "quasibinomial")

#marginal means 
library(emmeans)
B73_TIR.rg <- ref_grid(B73_TIR_glm)
B73_TIR_emmeans_loc <- emmeans(B73_TIR.rg, "genomic_loc", infer = T, level = 0.95)
pairs(B73_TIR_emmeans_loc)
B73_TIR_emmeans_sub <- emmeans(B73_TIR.rg, "subgenome", infer = T, level = 0.95)
pairs(B73_TIR_emmeans_sub)
B73_TIR_emmeans_nest <- emmeans(B73_TIR.rg, "nested_status", infer = T, level = 0.95)
pairs(B73_TIR_emmeans_nest)
#B73 TIRs only non-nested TEs, for modeling methylation data
B73_TIRs_nonnested <-subset(B73_TIRs, nested_status == "non_nested_TE")
B73_TIR_nonnested_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_TIRs_nonnested, family = "quasibinomial")
B73_TIR_nonnested_indvl_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high + chh_high + chg_high, data = B73_TIRs_nonnested, family = "quasibinomial")

#just B73 Helitrons
B73_Helitrons <- subset(WiDiv508_refB73_TEmetadata_low.ambig, order == "Helitron")
B73_Helitrons_glm <- glm(prop_present ~ factor(genomic_loc) + factor(subgenome) + nested_status + TE_len + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_Helitrons, family  = "quasibinomial")
B73_Helitron.rg <- ref_grid(B73_Helitrons_glm)
B73_Helitron_emmeans_loc <- emmeans(B73_Helitron.rg, "genomic_loc", infer = T, level = 0.95)
pairs(B73_Helitron_emmeans_loc)
B73_Helitron_emmeans_sub <- emmeans(B73_Helitron.rg, "subgenome", infer = T, level = 0.95)
pairs(B73_Helitron_emmeans_sub)
B73_Helitron_emmeans_nest <- emmeans(B73_Helitron.rg, "nested_status", infer = T, level = 0.95)
pairs(B73_Helitron_emmeans_nest)
#B73 Helitrons only non-nested TEs, for modeling methylation data
B73_Helitrons_nonnested <-subset(B73_Helitrons, nested_status == "non_nested_TE")
B73_Helitrons_nonnested_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_Helitrons_nonnested, family = "quasibinomial")
B73_Helitrons_nonnested_indvl_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high + chh_high + chg_high, data = B73_Helitrons_nonnested, family = "quasibinomial")

#Just B73 LTRS
B73_LTRs_glm <- glm(prop_present ~ factor(genomic_loc) + factor(subgenome) + nested_status + TE_len + LTR_similarity + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = WiDiv508_refB73_TEmetadata_LTRs_low.ambig, family = "quasibinomial")
summary(B73_LTRs_glm)
B73_LTR.rg <- ref_grid(B73_LTRs_glm)
B73_LTR_emmeans_loc <- emmeans(B73_LTR.rg, "genomic_loc", infer = T, level = 0.95)
pairs(B73_LTR_emmeans_loc)
B73_LTR_emmeans_sub <- emmeans(B73_LTR.rg, "subgenome", infer = T, level = 0.95)
pairs(B73_LTR_emmeans_sub)
B73_LTR_emmeans_nest <- emmeans(B73_LTR.rg, "nested_status", infer = T, level = 0.95)
pairs(B73_LTR_emmeans_nest)
#B73 LTRS only non-nested TEs, for modeling methylation data
B73_LTRs_nonnested <-subset(WiDiv508_refB73_TEmetadata_LTRs_low.ambig, nested_status == "non_nested_TE")
B73_LTRs_nonnested_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_LTRs_nonnested, family = "quasibinomial")
B73_LTRs_nonnested_indvl_glm <- glm(prop_present ~ genomic_loc + subgenome + TE_len + fam_size + cg_high + chh_high + chg_high, data = B73_LTRs_nonnested, family = "quasibinomial")


#glm for B73 TEs witin 10 kb of a gene
#genomic_loc and location_up_neg have the same info in different forms 
#fix bins_2kb errors 
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb == "TE_1kb_upstream"] <- "up_0kb.2kb"
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb == "TE_1.5kb_downstream"] <- "down_0kb.2kb"
WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb[WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig$bins_2kb == "TE_1kb_downstream"] <- "down_0kb.2kb"

B73_within10kb_TIR <- subset(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, order = "TIR")
B73_within10kb_Helitron <- subset(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, order = "Helitron")
B73_within10kb_LTR <- subset(WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, order = "LTR")
#add LTR similarity data 
WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig_sim <- data.frame(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig$TE_name, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig$LTR_similarity)
names(WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig_sim) <- c("TE_name", "LTR_similarity")
B73_within10kb_LTR <- merge(B73_within10kb_LTR, WiDiv508_refB73_TEmetadata_10kb_genes_LTRs_low.ambig_sim, by = "TE_name")

B73_within10kb_glm_TIR <- glm(prop_present ~ subgenome + TE_len + as.numeric(distance_to_gene) + fam_size + exp_prop_gene + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_within10kb_TIR, family = "quasibinomial")
summary(B73_within10kb_glm_TIR)

#B73_within10kb_glm_TIR_locfactor <- glm(prop_present ~ subgenome + bins_2kb + TE_len + fam_size + exp_prop_gene + cg_high_fam + chh_high_fam + chg_high_fam, data = B73_within10kb_TIR, family = "quasibinomial")
#summary(B73_within10kb_glm_TIR_locfactor)

B73_within10kb_glm_Helitron <- glm(prop_present ~ subgenome + TE_len + as.numeric(distance_to_gene) + exp_prop_gene + fam_size + exp_prop_gene + cg_high + chh_high + chg_high, data = B73_within10kb_Helitron, family = "quasibinomial")
summary(B73_within10kb_glm_Helitron)

#B73_within10kb_glm_Helitron_locfactor <- glm(prop_present ~ subgenome + bins_2kb + TE_len + exp_prop_gene + fam_size + exp_prop_gene + cg_high + chh_high + chg_high, data = B73_within10kb_Helitron, family = "quasibinomial")
#summary(B73_within10kb_glm_Helitron_locfactor)

#LTRS
B73_within10kb_glm_LTRs <- glm(prop_present ~ subgenome + as.numeric(distance_to_gene) + TE_len + exp_prop_gene + LTR_similarity + fam_size + cg_high + chh_high + chg_high, data = B73_within10kb_LTR, family = "quasibinomial")
summary(B73_within10kb_glm_LTRs)

B73_within10kb_glm_LTRs_locfactor <- glm(prop_present ~ subgenome + bins_2kb + TE_len + exp_prop_gene + LTR_similarity + fam_size + cg_high + chh_high + chg_high, data = B73_within10kb_LTR, family = "quasibinomial")
summary(B73_within10kb_glm_LTRs_locfactor)


#write out files with correctly formatted data 
write.table(x = WiDiv508_refB73_TEmetadata_low.ambig, file = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEmetadata_allmodeldata.txt", sep = "\t", row.names = F)
write.table(x = WiDiv508_refB73_TEmetadata_10kb_genes_low.ambig, file = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_TEmetadata_10kbgenes_allmodeldata.txt", sep = "\t", row.names = F)

write.table(x = WiDiv508_refB73_TEmetadata_LTRs_low.ambig, file = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_LTR_TEmetadata_allmodeldata.txt", sep = "\t", row.names = F)
write.table(x = B73_within10kb_LTR, file = "~/Documents/TE_project/ref_B73_TEcalling/WiDiv508_analysis/WiDiv508_refB73_LTR_TEmetadata_10kbgenes_allmodeldata.txt", sep = "\t", row.names = F)


#end of script 

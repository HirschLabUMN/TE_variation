#format data and make graphs to illustrate TE calling rates using bedtools 
library(ggplot2)
library(reshape2)
library(tidyverse)
library(plyr)


#small function to call true/false calling rates; specific to this data frame format
calc_true_false_pos <- function(dataframe, genotype, bp_len) {
	#dataframe is output of calc_sensitivity_cov_counts function
	#this is specific to the above files
	cov_lim <- c(2,3,4,5,6,7,8)
	empty_df <- data.frame(population = character(), cov_limit = numeric(), bp_length = numeric(), total_present = numeric(), total_absent = numeric(), true_present_calls = numeric(), true_absent_calls = numeric(), false_pres_calls = numeric(), false_absent_calls = numeric(), true_pres_rate_noamb = numeric(), true_abs_rate_noamb = numeric(), false_pres_rate_noamb = numeric(), false_abs_rate_noamb = numeric(), n_ambiguous = numeric())
	for (num in cov_lim) {
		all_present <- nrow(dataframe[which(dataframe$TE_status == "Present"),])
		all_absent <- nrow(dataframe[which(dataframe$TE_status == "Absent"),]) + nrow(dataframe[which(dataframe$TE_status == "Absent.SD"),])
		current_col = num + 5 #get data from the right column
		#numbers specific to each coverage limit
		present_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),])
		true_present <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),])
		absent_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])
		true_absent <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),])
		true_present_rate <- true_present/present_tes #actually called present, over all present, non-ambiguous calls
		true_absent_rate <- true_absent/absent_tes #actually called absent, over all absent, non-ambiguous calls
		false_absent_rate <- nrow(dataframe[which(dataframe[,current_col]== "call_false_absent"),])/present_tes #false absent = called as absent, but actually present
		false_present_rate <- nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])/absent_tes #false present = called as present, but actually absent
		n_ambiguous <- nrow(dataframe[which(dataframe[,current_col] == "ambiguous"),]) 
		empty_df <- rbind(empty_df, data.frame(population = genotype, cov_limit = num, bp_length = bp_len, total_present = all_present, total_absent = all_absent, true_present_calls = true_present, true_absent_calls = true_absent, false_pres_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),]), false_absent_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),]), true_pres_rate_noamb = true_present_rate, true_abs_rate_noamb = true_absent_rate, false_pres_rate_noamb = false_present_rate, false_abs_rate_noamb = false_absent_rate, n_ambiguous = n_ambiguous))
	}	
	return(empty_df)
}

calc_true_false_pos_lowcov <- function(dataframe, genotype, bp_len) {
	#dataframe is output of calc_sensitivity_cov_counts function
	#this is specific to the above files
	cov_lim <- c(1,2,3,4,5,6,7,8)
	empty_df <- data.frame(population = character(), cov_limit = numeric(), bp_length = numeric(), total_present = numeric(), total_absent = numeric(), true_present_calls = numeric(), true_absent_calls = numeric(), false_pres_calls = numeric(), false_absent_calls = numeric(), true_pres_rate_noamb = numeric(), true_abs_rate_noamb = numeric(), false_pres_rate_noamb = numeric(), false_abs_rate_noamb = numeric(), n_ambiguous = numeric())
	for (num in cov_lim) {
		all_present <- nrow(dataframe[which(dataframe$TE_status == "Present"),])
		all_absent <- nrow(dataframe[which(dataframe$TE_status == "Absent"),]) + nrow(dataframe[which(dataframe$TE_status == "Absent.SD"),])
		current_col = num + 6 #get data from the right column
		#numbers specific to each coverage limit
		present_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),])
		true_present <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),])
		absent_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])
		true_absent <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),])
		true_present_rate <- true_present/present_tes #actually called present, over all present, non-ambiguous calls
		true_absent_rate <- true_absent/absent_tes #actually called absent, over all absent, non-ambiguous calls
		false_absent_rate <- nrow(dataframe[which(dataframe[,current_col]== "call_false_absent"),])/present_tes #false absent = called as absent, but actually present
		false_present_rate <- nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])/absent_tes #false present = called as present, but actually absent
		n_ambiguous <- nrow(dataframe[which(dataframe[,current_col] == "ambiguous"),]) 
		empty_df <- rbind(empty_df, data.frame(population = genotype, cov_limit = num, bp_length = bp_len, total_present = all_present, total_absent = all_absent, true_present_calls = true_present, true_absent_calls = true_absent, false_pres_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),]), false_absent_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),]), true_pres_rate_noamb = true_present_rate, true_abs_rate_noamb = true_absent_rate, false_pres_rate_noamb = false_present_rate, false_abs_rate_noamb = false_absent_rate, n_ambiguous = n_ambiguous))
	}	
	return(empty_df)
}

#re-do true true/false positive calculations for Mo17 10 bp subsample .45
Mo17_calls_10bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.45_TEcov_10bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.45_calls <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.45, "Mo17", 10)

calc_true_false_ratio <- function(dataframe, genotype, bp_len) {
	#dataframe is output of calc_sensitivity_cov_counts function
	#this is specific to the above files
	cov_lim <- c(0.25,0.40,0.5,0.7,0.8)
	count = 1
	empty_df <- data.frame(population = character(), cov_limit = numeric(), bp_length = numeric(), total_present = numeric(), total_absent = numeric(), true_present_calls = numeric(), true_absent_calls = numeric(), false_pres_calls = numeric(), false_absent_calls = numeric(), true_pres_rate_noamb = numeric(), true_abs_rate_noamb = numeric(), false_pres_rate_noamb = numeric(), false_abs_rate_noamb = numeric(), n_ambiguous = numeric())
	
	for (num in cov_lim) {
		all_present <- nrow(dataframe[which(dataframe$TE_status == "Present"),])
		all_absent <- nrow(dataframe[which(dataframe$TE_status == "Absent"),]) + nrow(dataframe[which(dataframe$TE_status == "Absent.SD"),])
		current_col = count + 8 #get data from the right column
		#numbers specific to each coverage limit
		present_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),])
		true_present <- nrow(dataframe[which(dataframe[,current_col] == "call_present"),])
		absent_tes <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),]) + nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])
		true_absent <- nrow(dataframe[which(dataframe[,current_col] == "call_absent"),])
		true_present_rate <- true_present/present_tes #actually called present, over all present, non-ambiguous calls
		true_absent_rate <- true_absent/absent_tes #actually called absent, over all absent, non-ambiguous calls
		false_absent_rate <- nrow(dataframe[which(dataframe[,current_col]== "call_false_absent"),])/present_tes #false absent = called as absent, but actually present
		false_present_rate <- nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),])/absent_tes #false present = called as present, but actually absent
		n_ambiguous <- nrow(dataframe[which(dataframe[,current_col] == "ambiguous"),]) 
		empty_df <- rbind(empty_df, data.frame(population = genotype, cov_limit = num, bp_length = bp_len, total_present = all_present, total_absent = all_absent, true_present_calls = true_present, true_absent_calls = true_absent, false_pres_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_present"),]), false_absent_calls = nrow(dataframe[which(dataframe[,current_col] == "call_false_absent"),]), true_pres_rate_noamb = true_present_rate, true_abs_rate_noamb = true_absent_rate, false_pres_rate_noamb = false_present_rate, false_abs_rate_noamb = false_absent_rate, n_ambiguous = n_ambiguous))
		current_col = current_col + 1
	}	
	return(empty_df)
}

#read in files that have the coverage calls 
Mo17_calls_1bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_1bp.tf.pos.txt")
PH207_calls_1bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_1bp.tf.pos.txt")
W22_calls_1bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_1bp.tf.pos.txt")
Mo17_calls_3bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_3bp.tf.pos.txt")
PH207_calls_3bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_3bp.tf.pos.txt")
W22_calls_3bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_3bp.tf.pos.txt")
Mo17_calls_5bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_5bp.tf.pos.txt")
PH207_calls_5bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_5bp.tf.pos.txt")
W22_calls_5bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_5bp.tf.pos.txt")
Mo17_calls_10bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov.tf.pos.txt")
PH207_calls_10bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov.tf.pos.txt")
W22_calls_10bp <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov.tf.pos.txt")
#coverage rates using inner coverage values 
Mo17_calls_1bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_1bp.tf.pos_inner.cov.txt")
PH207_calls_1bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_1bp.tf.pos_inner.cov.txt")
W22_calls_1bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_1bp.tf.pos_inner.cov.txt")
Mo17_calls_3bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_3bp.tf.pos_inner.cov.txt")
PH207_calls_3bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_3bp.tf.pos_inner.cov.txt")
W22_calls_3bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_3bp.tf.pos_inner.cov.txt")
Mo17_calls_5bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_5bp.tf.pos.txt")
PH207_calls_5bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_5bp.tf.pos_inner.cov.txt")
W22_calls_5bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_5bp.tf.pos_inner.cov.txt")
Mo17_calls_10bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_10bp.tf.pos_inner.cov.txt")
PH207_calls_10bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_10bp.tf.pos_inner.cov.txt")
W22_calls_10bp_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_10bp.tf.pos_inner.cov.txt")
#read in subsampled bams with correct 1 bp calls 
Mo17_calls_1bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.20_TEcov_1bp.tf.pos.inner.txt")
Mo17_calls_3bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.20_TEcov_3bp.tf.pos.inner.txt")
Mo17_calls_5bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.20_TEcov_5bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.20_TEcov_10bp.tf.pos.inner.txt")
PH207_calls_1bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.20_TEcov_1bp.tf.pos.inner.txt")
PH207_calls_3bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.20_TEcov_3bp.tf.pos.inner.txt")
PH207_calls_5bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.20_TEcov_5bp.tf.pos.innertxt")
PH207_calls_10bp_inner_sub0.20 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.20_TEcov_10bp.tf.pos.inner.txt")
Mo17_calls_1bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.32_TEcov_1bp.tf.pos.inner.txt")
Mo17_calls_3bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.32_TEcov_3bp.tf.pos.inner.txt")
Mo17_calls_5bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.32_TEcov_5bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.32_TEcov_10bp.tf.pos.inner.txt")
PH207_calls_1bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.32_TEcov_1bp.tf.pos.inner.txt")
PH207_calls_3bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.32_TEcov_3bp.tf.pos.inner.txt")
PH207_calls_5bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.32_TEcov_5bp.tf.pos.inner.txt")
PH207_calls_10bp_inner_sub0.32 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.32_TEcov_10bp.tf.pos.inner.txt")
Mo17_calls_1bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.45_TEcov_1bp.tf.pos.inner.txt")
Mo17_calls_3bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.45_TEcov_3bp.tf.pos.inner.txt")
Mo17_calls_5bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.45_TEcov_5bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.45_TEcov_10bp.tf.pos.inner.txt")
PH207_calls_1bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.45_TEcov_1bp.tf.pos.inner.txt")
PH207_calls_3bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.45_TEcov_3bp.tf.pos.inner.txt")
PH207_calls_5bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.45_TEcov_5bp.tf.pos.inner.txt")
PH207_calls_10bp_inner_sub0.45 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.45_TEcov_10bp.tf.pos.inner.txt")
W22_calls_1bp_1bp_accurate <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcov_1bp.tf.pos.inner_1pb.lim_accurate.txt")
W22_calls_3bp_1bp_accurate <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcov_3bp.tf.pos.inner_1pb.lim_accurate.txt")
W22_calls_5bp_1bp_accurate <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcov_5bp.tf.pos.inner_1pb.lim_accurate.txt")
W22_calls_10bp_1bp_accurate <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcov_10bp.tf.pos.inner_1pb.lim_accurate.txt")
#W22 clipped reads 
#W22_calls_10bp_clipped <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_clipped.reads_TEcov_10bp.tf.pos.inner.txt")
#read in 2 category reads, 4 sides 
Mo17_call_1bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_1bp.tf.pos_2cat.txt")
Mo17_call_3bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_3bp.tf.pos_2cat.txt")
Mo17_call_5bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_5bp.tf.pos_2cat.txt")
Mo17_call_10bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_10bp.tf.pos_2cat.txt")
PH207_call_1bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_1bp.tf.pos_2cat.txt")
PH207_call_3bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_3bp.tf.pos_2cat.txt")
PH207_call_5bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_5bp.tf.pos_2cat.txt")
PH207_call_10bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_10bp.tf.pos_2cat.txt")
W22_call_1bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_1bp.tf.pos_2cat.txt")
W22_call_3bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_3bp.tf.pos_2cat.txt")
W22_call_5bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_5bp.tf.pos_2cat.txt")
W22_call_10bp_2cat_4sides <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_10bp.tf.pos_2cat.txt")
#read in 2 category reads, inner sides of TE ends 
Mo17_call_1bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_1bp.tf.pos_2cat.inner.cov.txt")
Mo17_call_3bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_3bp.tf.pos_2cat.inner.cov.txt")
Mo17_call_5bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_5bp.tf.pos_2cat.inner.cov.txt")
Mo17_call_10bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_10bp.tf.pos_2cat.inner.cov.txt")
PH207_call_1bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_1bp.tf.pos_2cat.inner.cov.txt")
PH207_call_3bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_3bp.tf.pos_2cat.inner.cov.txt")
PH207_call_5bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_5bp.tf.pos_2cat.inner.cov.txt")
PH207_call_10bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_10bp.tf.pos_2cat.inner.cov.txt")
W22_call_1bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_1bp.tf.pos_2cat.inner.cov.txt")
W22_call_3bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_3bp.tf.pos_2cat.inner.cov.txt")
W22_call_5bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_5bp.tf.pos_2cat.inner.cov.txt")
W22_call_10bp_2cat_inner <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_10bp.tf.pos_2cat.inner.cov.txt")
#more subsample files 
Mo17_calls_1bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.66_TEcov_1bp.tf.pos.inner.txt")
Mo17_calls_3bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.66_TEcov_3bp.tf.pos.inner.txt")
Mo17_calls_5bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.66_TEcov_5bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.66_TEcov_10bp.tf.pos.inner.txt")
PH207_calls_1bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.66_TEcov_1bp.tf.pos.inner.txt")
PH207_calls_3bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.66_TEcov_3bp.tf.pos.inner.txt")
PH207_calls_5bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.66_TEcov_5bp.tf.pos.inner.txt")
PH207_calls_10bp_inner_sub0.66 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.66_TEcov_10bp.tf.pos.inner.txt")
Mo17_calls_1bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.83_TEcov_1bp.tf.pos.inner.txt")
Mo17_calls_3bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.83_TEcov_3bp.tf.pos.inner.txt")
Mo17_calls_5bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.83_TEcov_5bp.tf.pos.inner.txt")
Mo17_calls_10bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73v4_subsample0.83_TEcov_10bp.tf.pos.inner.txt")
PH207_calls_1bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.83_TEcov_1bp.tf.pos.inner.txt")
PH207_calls_3bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.83_TEcov_3bp.tf.pos.inner.txt")
PH207_calls_5bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.83_TEcov_5bp.tf.pos.inner.txt")
PH207_calls_10bp_inner_sub0.83 <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73v4_subsample0.83_TEcov_10bp.tf.pos.inner.txt")


#edit subsample calling rates 
#Why? 1 bp calling was messed up, so I'm fixing that here 
cov1_calling_rates <- function(subsample_df) {
	#correct calls_1bp_lim column entries
	for (row in 1:nrow(subsample_df)){
		if(subsample_df[row,2] == "Present" & subsample_df[row,4] >= 1 & subsample_df[row,5] >= 1) {
			subsample_df[row,7] = "call_present"
		} else if(subsample_df[row,2] == "Absent" & subsample_df[row,4] == 0 & subsample_df[row,5] == 0 | subsample_df[row,2] == "Absent.SD" & subsample_df[row,4] == 0 & subsample_df[row,5] == 0) {
			subsample_df[row,7] = "call_absent"
		} else if (subsample_df[row,2] == "Present" & subsample_df[row,4] == 0 & subsample_df[row,5] == 0) {
			subsample_df[row,7] = "call_false_absent"
		} else if (subsample_df[row,2] == "Absent" & subsample_df[row,4] >= 1 & subsample_df[row,5] >= 1 | subsample_df[row,2] == "Absent.SD" & subsample_df[row,4] >= 1 & subsample_df[row,5] >= 1) {
			subsample_df[row,7] = "call_false_present"
		} else if (subsample_df[row,2] == "Unresolved") {
			subsample_df[row,7] = "undefined"
		} else {
			subsample_df[row,7] = "ambiguous"
		}
	}
	return(subsample_df) 
} #this was run on mesabi
#ratio files 
#Mo17_calls_1bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_1bp.tf.pos_ratio.txt")
#Mo17_calls_3bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_3bp.tf.pos_ratio.txt")
#Mo17_calls_5bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_5bp.tf.pos_ratio.txt")
#Mo17_calls_10bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/Mo17_refB73_TEcov_10bp.tf.pos_ratio.txt")
#PH207_calls_1bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_1bp.tf.pos_ratio.txt")
#PH207_calls_3bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_3bp.tf.pos_ratio.txt")
#PH207_calls_5bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_5bp.tf.pos_ratio.txt")
#PH207_calls_10bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/PH207_refB73_TEcov_10bp.tf.pos_ratio.txt")
#W22_calls_1bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_1bp.tf.pos_ratio.txt")
#W22_calls_3bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_3bp.tf.pos_ratio.txt")
#W22_calls_5bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_5bp.tf.pos_ratio.txt")
#W22_calls_10bp_ratio <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73_TEcov_10bp.tf.pos_ratio.txt")

#get accuracy rates 
#4 points, full coverage
Mo17_1bp_calling_rates <- calc_true_false_pos(Mo17_calls_1bp, "Mo17", 1)
PH207_1bp_calling_rates <- calc_true_false_pos(PH207_calls_1bp, "PH207", 1)
W22_1bp_calling_rates <- calc_true_false_pos(W22_calls_1bp, "W22", 1)
Mo17_3bp_calling_rates <- calc_true_false_pos(Mo17_calls_3bp, "Mo17", 3)
PH207_3bp_calling_rates <- calc_true_false_pos(PH207_calls_3bp, "PH207", 3)
W22_3bp_calling_rates <- calc_true_false_pos(W22_calls_3bp, "W22", 3)
Mo17_5bp_calling_rates <- calc_true_false_pos(Mo17_calls_5bp, "Mo17", 5)
PH207_5bp_calling_rates <- calc_true_false_pos(PH207_calls_5bp, "PH207", 5)
W22_5bp_calling_rates <- calc_true_false_pos(W22_calls_5bp, "W22", 5)
Mo17_10bp_calling_rates <- calc_true_false_pos(Mo17_calls_10bp, "Mo17", 10)
PH207_10bp_calling_rates <- calc_true_false_pos(PH207_calls_10bp, "PH207", 10)
W22_10bp_calling_rates <- calc_true_false_pos(W22_calls_10bp, "W22", 10)
TEcalling_rates_4points <- do.call("rbind", list(Mo17_1bp_calling_rates, PH207_1bp_calling_rates, W22_1bp_calling_rates, Mo17_3bp_calling_rates, PH207_3bp_calling_rates, W22_3bp_calling_rates, Mo17_5bp_calling_rates, PH207_5bp_calling_rates, W22_5bp_calling_rates, Mo17_10bp_calling_rates, PH207_10bp_calling_rates, W22_10bp_calling_rates))
write.table(TEcalling_rates_4points, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_4points.txt", sep = "\t", row.names = F)
#inner, full coverage
Mo17_1bp_calling_rates_inner <- calc_true_false_pos(Mo17_calls_1bp_inner, "Mo17", 1)
PH207_1bp_calling_rates_inner <- calc_true_false_pos(PH207_calls_1bp_inner, "PH207", 1)
W22_1bp_calling_rates_inner <- calc_true_false_pos(W22_calls_1bp_inner, "W22", 1)
Mo17_3bp_calling_rates_inner <- calc_true_false_pos(Mo17_calls_3bp_inner, "Mo17", 3)
PH207_3bp_calling_rates_inner <- calc_true_false_pos(PH207_calls_3bp_inner, "PH207", 3)
W22_3bp_calling_rates_inner <- calc_true_false_pos(W22_calls_3bp_inner, "W22", 3)
Mo17_5bp_calling_rates_inner <- calc_true_false_pos(Mo17_calls_5bp_inner, "Mo17", 5)
PH207_5bp_calling_rates_inner <- calc_true_false_pos(PH207_calls_5bp_inner, "PH207", 5)
W22_5bp_calling_rates_inner <- calc_true_false_pos(W22_calls_5bp_inner, "W22", 5)
Mo17_10bp_calling_rates_inner <- calc_true_false_pos(Mo17_calls_10bp_inner, "Mo17", 10)
PH207_10bp_calling_rates_inner <- calc_true_false_pos(PH207_calls_10bp_inner, "PH207", 10)
W22_10bp_calling_rates_inner <- calc_true_false_pos(W22_calls_10bp_inner, "W22", 10)
TEcalling_rates_inner <- do.call("rbind", list(Mo17_1bp_calling_rates_inner, PH207_1bp_calling_rates_inner, W22_1bp_calling_rates_inner, Mo17_3bp_calling_rates_inner, PH207_3bp_calling_rates_inner, W22_3bp_calling_rates_inner, Mo17_5bp_calling_rates_inner, PH207_5bp_calling_rates_inner, W22_5bp_calling_rates_inner, Mo17_10bp_calling_rates_inner, PH207_10bp_calling_rates_inner, W22_10bp_calling_rates_inner))
write.table(TEcalling_rates_inner, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner.txt", sep = "\t", row.names = F)
#summarize Mo17 and PH207 calls 
TEcalling_rates_Mo17PH207 <- do.call("rbind", list(Mo17_1bp_calling_rates, PH207_1bp_calling_rates, Mo17_3bp_calling_rates, PH207_3bp_calling_rates, Mo17_5bp_calling_rates, PH207_5bp_calling_rates, Mo17_10bp_calling_rates, PH207_10bp_calling_rates))
Mo17_PH207_summary <- ddply(TEcalling_rates_Mo17PH207, ~cov_limit + bp_length, summarise, true_pres_rate = mean(true_pres_rate_noamb), true_abs_rate = mean(true_abs_rate_noamb), nambig = mean(n_ambiguous))
write.table(Mo17_PH207_summary, "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_4points_Mo17PH207_mean.txt", sep = "\t", row.names = F)
TEcalling_rates_Mo17PH207_inner <- do.call("rbind", list(Mo17_1bp_calling_rates_inner, PH207_1bp_calling_rates_inner, Mo17_3bp_calling_rates_inner, PH207_3bp_calling_rates_inner, Mo17_5bp_calling_rates_inner, PH207_5bp_calling_rates_inner, Mo17_10bp_calling_rates_inner, PH207_10bp_calling_rates_inner))
Mo17_PH207_summary_inner <- ddply(TEcalling_rates_Mo17PH207_inner, ~ cov_limit + bp_length, summarise, true_pres_rate = mean(true_pres_rate_noamb), true_abs_rate = mean(true_abs_rate_noamb), nambig = mean(n_ambiguous))
write.table(Mo17_PH207_summary_inner, "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_Mo17PH207_mean.txt", sep = "\t", row.names = F)
#TE calling using raio method
#Mo17_1bp_calling_rates_ratio <- calc_true_false_ratio(Mo17_calls_1bp_ratio, "Mo17", 1)
#Mo17_3bp_calling_rates_ratio <- calc_true_false_ratio(Mo17_calls_3bp_ratio, "Mo17", 3)
#Mo17_5bp_calling_rates_ratio <- calc_true_false_ratio(Mo17_calls_5bp_ratio, "Mo17", 5)
#Mo17_10bp_calling_rates_ratio <- calc_true_false_ratio(Mo17_calls_10bp_ratio, "Mo17", 10)
#PH207_1bp_calling_rates_ratio <- calc_true_false_ratio(PH207_calls_1bp_ratio, "PH207", 1)
#PH207_3bp_calling_rates_ratio <- calc_true_false_ratio(PH207_calls_3bp_ratio, "PH207", 3)
#PH207_5bp_calling_rates_ratio <- calc_true_false_ratio(PH207_calls_5bp_ratio, "PH207", 5)
#PH207_10bp_calling_rates_ratio <- calc_true_false_ratio(PH207_calls_10bp_ratio, "PH207", 10)
#W22_1bp_calling_rates_ratio <- calc_true_false_ratio(W22_calls_1bp_ratio, "W22", 1)
#W22_3bp_calling_rates_ratio <- calc_true_false_ratio(W22_calls_3bp_ratio, "W22", 3)
#W22_5bp_calling_rates_ratio <- calc_true_false_ratio(W22_calls_5bp_ratio, "W22", 5)
#W22_10bp_calling_rates_ratio <- calc_true_false_ratio(W22_calls_10bp_ratio, "W22", 10)
#TEcalling_rates_ratio <- do.call("rbind", list(Mo17_1bp_calling_rates_ratio, Mo17_3bp_calling_rates_ratio, Mo17_5bp_calling_rates_ratio, Mo17_10bp_calling_rates_ratio, PH207_1bp_calling_rates_ratio, PH207_3bp_calling_rates_ratio, PH207_5bp_calling_rates_ratio, PH207_10bp_calling_rates_ratio, W22_1bp_calling_rates_ratio, W22_3bp_calling_rates_ratio, W22_5bp_calling_rates_ratio, W22_10bp_calling_rates_ratio))
#write.table(TEcalling_rates_ratio, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_ratio.txt", sep = "\t", row.names = F)
#TE calling using subsampled TEs, inner coverage
Mo17_1bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(Mo17_calls_1bp_inner_sub0.20, "Mo17", 1)
PH207_1bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(PH207_calls_1bp_inner_sub0.20, "PH207", 1)
Mo17_3bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(Mo17_calls_3bp_inner_sub0.20, "Mo17", 3)
PH207_3bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(PH207_calls_3bp_inner_sub0.20, "PH207", 3)
Mo17_5bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(Mo17_calls_5bp_inner_sub0.20, "Mo17", 5)
PH207_5bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(PH207_calls_5bp_inner_sub0.20, "PH207", 5)
Mo17_10bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.20, "Mo17", 10)
PH207_10bp_calling_rates_inner_sub0.20 <- calc_true_false_pos_lowcov(PH207_calls_10bp_inner_sub0.20, "PH207", 10)
TEcalling_rates_inner_sub0.20 <- do.call("rbind", list(Mo17_1bp_calling_rates_inner_sub0.20, PH207_1bp_calling_rates_inner_sub0.20, Mo17_3bp_calling_rates_inner_sub0.20, PH207_3bp_calling_rates_inner_sub0.20, Mo17_5bp_calling_rates_inner_sub0.20, PH207_5bp_calling_rates_inner_sub0.20, Mo17_10bp_calling_rates_inner_sub0.20, PH207_10bp_calling_rates_inner_sub0.20))
write.table(TEcalling_rates_inner_sub0.20, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.20.txt", sep = "\t", row.names = F)
Mo17_1bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(Mo17_calls_1bp_inner_sub0.32, "Mo17", 1)
PH207_1bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(PH207_calls_1bp_inner_sub0.32, "PH207", 1)
Mo17_3bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(Mo17_calls_3bp_inner_sub0.32, "Mo17", 3)
PH207_3bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(PH207_calls_3bp_inner_sub0.32, "PH207", 3)
Mo17_5bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(Mo17_calls_5bp_inner_sub0.32, "Mo17", 5)
PH207_5bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(PH207_calls_5bp_inner_sub0.32, "PH207", 5)
Mo17_10bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.32, "Mo17", 10)
PH207_10bp_calling_rates_inner_sub0.32 <- calc_true_false_pos_lowcov(PH207_calls_10bp_inner_sub0.32, "PH207", 10)
TEcalling_rates_inner_sub0.32 <- do.call("rbind", list(Mo17_1bp_calling_rates_inner_sub0.32, PH207_1bp_calling_rates_inner_sub0.32, Mo17_3bp_calling_rates_inner_sub0.32, PH207_3bp_calling_rates_inner_sub0.32, Mo17_5bp_calling_rates_inner_sub0.32, PH207_5bp_calling_rates_inner_sub0.32, Mo17_10bp_calling_rates_inner_sub0.32, PH207_10bp_calling_rates_inner_sub0.32))
write.table(TEcalling_rates_inner_sub0.32, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.32.txt", sep = "\t", row.names = F)
Mo17_1bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(Mo17_calls_1bp_inner_sub0.45, "Mo17", 1)
PH207_1bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(PH207_calls_1bp_inner_sub0.45, "PH207", 1)
Mo17_3bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(Mo17_calls_3bp_inner_sub0.45, "Mo17", 3)
PH207_3bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(PH207_calls_3bp_inner_sub0.45, "PH207", 3)
Mo17_5bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(Mo17_calls_5bp_inner_sub0.45, "Mo17", 5)
PH207_5bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(PH207_calls_5bp_inner_sub0.45, "PH207", 5)
Mo17_10bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.45, "Mo17", 10)
PH207_10bp_calling_rates_inner_sub0.45 <- calc_true_false_pos_lowcov(PH207_calls_10bp_inner_sub0.45, "PH207", 10)
TEcalling_rates_inner_sub0.45 <- do.call("rbind", list(Mo17_1bp_calling_rates_inner_sub0.45, PH207_1bp_calling_rates_inner_sub0.45, Mo17_3bp_calling_rates_inner_sub0.45, PH207_3bp_calling_rates_inner_sub0.45, Mo17_5bp_calling_rates_inner_sub0.45, PH207_5bp_calling_rates_inner_sub0.32, Mo17_10bp_calling_rates_inner_sub0.45, PH207_10bp_calling_rates_inner_sub0.45))
write.table(TEcalling_rates_inner_sub0.45, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.45.txt", sep = "\t", row.names = F)
#get accuracy rates for W22 1bp lim calls 
W22_1bp_calling_rates_1bplim <- calc_true_false_pos_lowcov(W22_calls_1bp_1bp_accurate, "W22", 1)
W22_3bp_calling_rates_1bplim <- calc_true_false_pos_lowcov(W22_calls_3bp_1bp_accurate, "W22", 3)
W22_5bp_calling_rates_1bplim <- calc_true_false_pos_lowcov(W22_calls_5bp_1bp_accurate, "W22", 5)
W22_10bp_calling_rates_1bplim <- calc_true_false_pos_lowcov(W22_calls_10bp_1bp_accurate, "W22", 10)
TEcalling_rates_W22_1bplim <- do.call("rbind", list(W22_1bp_calling_rates_1bplim, W22_3bp_calling_rates_1bplim, W22_5bp_calling_rates_1bplim, W22_10bp_calling_rates_1bplim))
write.table(TEcalling_rates_W22_1bplim, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_W22_10bp_1bp.lim.accurate.txt", sep = "\t", row.names = F)
#compare W22 clipped reads with normal reads 
#10 bp TE fragment length 
W22_10bp_clipped_calling_rates_inner <- calc_true_false_pos(W22_calls_10bp_clipped, "W22", 10)
#2 bp lim didn't work 
#re-do TE calling 
W22_calls_10bp_clipped_callrates <- calc_sensitivity_cov_counts_innercov(W22_calls_10bp_clipped)
write.table(W22_calls_10bp_clipped_callrates, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_clipped.reads_TEcov_10bp.tf.pos.inner_corrected.txt")
W22_10bp_clipped_calling_rates_inner <- calc_true_false_pos(W22_calls_10bp_clipped_callrates, "W22", 10)
W22_10bp_clipped_calling_rates_inner$read_len <- "clipped_100bp"
W22_10bp_calling_rates_inner <- calc_true_false_pos(W22_calls_10bp_inner, "W22", 10)
W22_10bp_calling_rates_inner$read_len <- "orig_150bp"

W22_10bp_calling_inner_normal.clipped <- do.call("rbind", list(W22_10bp_clipped_calling_rates_inner, W22_10bp_calling_rates_inner))
write.table(W22_10bp_calling_inner_normal.clipped, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcalling_rates_10bp_normal.clipped.txt", sep = "\t", row.names = F)

#2 category, 4 side accuracy rates 
Mo17_1bp_call_rates_2cat_4sides <- calc_true_false_pos(Mo17_call_1bp_2cat_4sides, "Mo17", 1)
Mo17_3bp_call_rates_2cat_4sides <- calc_true_false_pos(Mo17_call_3bp_2cat_4sides, "Mo17", 3)
Mo17_5bp_call_rates_2cat_4sides <- calc_true_false_pos(Mo17_call_5bp_2cat_4sides, "Mo17", 5)
Mo17_10bp_call_rates_2cat_4sides <- calc_true_false_pos(Mo17_call_10bp_2cat_4sides, "Mo17", 10)
PH207_1bp_call_rates_2cat_4sides <- calc_true_false_pos(PH207_call_1bp_2cat_4sides, "PH207", 1)
PH207_3bp_call_rates_2cat_4sides <- calc_true_false_pos(PH207_call_3bp_2cat_4sides, "PH207", 3)
PH207_5bp_call_rates_2cat_4sides <- calc_true_false_pos(PH207_call_5bp_2cat_4sides, "PH207", 5)
PH207_10bp_call_rates_2cat_4sides <- calc_true_false_pos(PH207_call_10bp_2cat_4sides, "PH207", 10)
W22_1bp_call_rates_2cat_4sides <- calc_true_false_pos(W22_call_1bp_2cat_4sides, "W22", 1)
W22_3bp_call_rates_2cat_4sides <- calc_true_false_pos(W22_call_3bp_2cat_4sides, "W22", 3)
W22_5bp_call_rates_2cat_4sides <- calc_true_false_pos(W22_call_5bp_2cat_4sides, "W22", 5)
W22_10bp_call_rates_2cat_4sides <- calc_true_false_pos(W22_call_10bp_2cat_4sides, "W22", 10)
TEcalling_rates_2cat_4points <- do.call("rbind", list(Mo17_1bp_call_rates_2cat_4sides, PH207_1bp_call_rates_2cat_4sides, W22_1bp_call_rates_2cat_4sides, Mo17_3bp_call_rates_2cat_4sides, PH207_3bp_call_rates_2cat_4sides, W22_3bp_call_rates_2cat_4sides, Mo17_5bp_call_rates_2cat_4sides, PH207_5bp_call_rates_2cat_4sides, W22_5bp_call_rates_2cat_4sides, Mo17_10bp_call_rates_2cat_4sides, PH207_10bp_call_rates_2cat_4sides, W22_10bp_call_rates_2cat_4sides))
write.table(TEcalling_rates_2cat_4points, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_4points.txt", sep = "\t", row.names = F)
#2 category, 2 side accuracy rates 
Mo17_1bp_call_rates_2cat_inner <- calc_true_false_pos(Mo17_call_1bp_2cat_inner, "Mo17", 1)
Mo17_3bp_call_rates_2cat_inner <- calc_true_false_pos(Mo17_call_3bp_2cat_inner, "Mo17", 3)
Mo17_5bp_call_rates_2cat_inner <- calc_true_false_pos(Mo17_call_5bp_2cat_inner, "Mo17", 5)
Mo17_10bp_call_rates_2cat_inner <- calc_true_false_pos(Mo17_call_10bp_2cat_inner, "Mo17", 10)
PH207_1bp_call_rates_2cat_inner <- calc_true_false_pos(PH207_call_1bp_2cat_inner, "PH207", 1)
PH207_3bp_call_rates_2cat_inner <- calc_true_false_pos(PH207_call_3bp_2cat_inner, "PH207", 3)
PH207_5bp_call_rates_2cat_inner <- calc_true_false_pos(PH207_call_5bp_2cat_inner, "PH207", 5)
PH207_10bp_call_rates_2cat_inner <- calc_true_false_pos(PH207_call_10bp_2cat_inner, "PH207", 10)
W22_1bp_call_rates_2cat_inner <- calc_true_false_pos(W22_call_1bp_2cat_inner, "W22", 1)
W22_3bp_call_rates_2cat_inner <- calc_true_false_pos(W22_call_3bp_2cat_inner, "W22", 3)
W22_5bp_call_rates_2cat_inner <- calc_true_false_pos(W22_call_5bp_2cat_inner, "W22", 5)
W22_10bp_call_rates_2cat_inner <- calc_true_false_pos(W22_call_10bp_2cat_inner, "W22", 10)
TEcalling_rates_2cat_inner <- do.call("rbind", list(Mo17_1bp_call_rates_2cat_inner, PH207_1bp_call_rates_2cat_inner, W22_1bp_call_rates_2cat_inner, Mo17_3bp_call_rates_2cat_inner, PH207_3bp_call_rates_2cat_inner, W22_3bp_call_rates_2cat_inner, Mo17_5bp_call_rates_2cat_inner, PH207_5bp_call_rates_2cat_inner, W22_5bp_call_rates_2cat_inner, Mo17_10bp_call_rates_2cat_inner, PH207_10bp_call_rates_2cat_inner, W22_10bp_call_rates_2cat_inner))
write.table(TEcalling_rates_2cat_inner, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_inner.txt", sep = "\t", row.names = F)
#TE calling rates more subsample
Mo17_1bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(Mo17_calls_1bp_inner_sub0.66, "Mo17", 1)
PH207_1bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(PH207_calls_1bp_inner_sub0.66, "PH207", 1)
Mo17_3bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(Mo17_calls_3bp_inner_sub0.66, "Mo17", 3)
PH207_3bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(PH207_calls_3bp_inner_sub0.66, "PH207", 3)
Mo17_5bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(Mo17_calls_5bp_inner_sub0.66, "Mo17", 5)
PH207_5bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(PH207_calls_5bp_inner_sub0.66, "PH207", 5)
Mo17_10bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.66, "Mo17", 10)
PH207_10bp_calling_rates_inner_sub0.66 <- calc_true_false_pos_lowcov(PH207_calls_10bp_inner_sub0.66, "PH207", 10)
TEcalling_rates_inner_sub0.66 <- do.call("rbind", list(Mo17_1bp_calling_rates_inner_sub0.66, PH207_1bp_calling_rates_inner_sub0.66, Mo17_3bp_calling_rates_inner_sub0.66, PH207_3bp_calling_rates_inner_sub0.66, Mo17_5bp_calling_rates_inner_sub0.66, PH207_5bp_calling_rates_inner_sub0.66, Mo17_10bp_calling_rates_inner_sub0.66, PH207_10bp_calling_rates_inner_sub0.66))
write.table(TEcalling_rates_inner_sub0.66, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.66.txt", sep = "\t", row.names = F)
Mo17_1bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(Mo17_calls_1bp_inner_sub0.83, "Mo17", 1)
PH207_1bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(PH207_calls_1bp_inner_sub0.83, "PH207", 1)
Mo17_3bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(Mo17_calls_3bp_inner_sub0.83, "Mo17", 3)
PH207_3bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(PH207_calls_3bp_inner_sub0.83, "PH207", 3)
Mo17_5bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(Mo17_calls_5bp_inner_sub0.83, "Mo17", 5)
PH207_5bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(PH207_calls_5bp_inner_sub0.83, "PH207", 5)
Mo17_10bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(Mo17_calls_10bp_inner_sub0.83, "Mo17", 10)
PH207_10bp_calling_rates_inner_sub0.83 <- calc_true_false_pos_lowcov(PH207_calls_10bp_inner_sub0.83, "PH207", 10)
TEcalling_rates_inner_sub0.83 <- do.call("rbind", list(Mo17_1bp_calling_rates_inner_sub0.83, PH207_1bp_calling_rates_inner_sub0.83, Mo17_3bp_calling_rates_inner_sub0.83, PH207_3bp_calling_rates_inner_sub0.83, Mo17_5bp_calling_rates_inner_sub0.83, PH207_5bp_calling_rates_inner_sub0.83, Mo17_10bp_calling_rates_inner_sub0.83, PH207_10bp_calling_rates_inner_sub0.83))
write.table(TEcalling_rates_inner_sub0.83, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.83.txt", sep = "\t", row.names = F)



#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#read in summary files to avoid re-doing all the above again and again

#4 points re-shape data frame
TEcalling_rates_4points <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_4points.txt")
TEcalling_rates_4points_covsample <- subset(TEcalling_rates_4points, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_4points_melt_rates <- melt(TEcalling_rates_4points, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_4points_melt_rates_covsample <- subset(TEcalling_rates_4points_melt_rates, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_4points_melt_numbers <- melt(TEcalling_rates_4points, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_present_calls", "true_absent_calls", "false_pres_calls", "false_absent_calls"))
TEcalling_rates_4points_melt_numbers_covsample <- subset(TEcalling_rates_4points_melt_numbers, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#inner coverage; re-shape data frame
TEcalling_rates_inner <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner.txt")
TEcalling_rates_inner_covsample <- subset(TEcalling_rates_inner, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_inner_melt_rates <- melt(TEcalling_rates_inner, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_inner_melt_rates_covsample <- subset(TEcalling_rates_inner_melt_rates, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_inner_melt_numbers <- melt(TEcalling_rates_inner, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_present_calls", "true_absent_calls", "false_pres_calls", "false_absent_calls"))
TEcalling_rates_inner_melt_numbers_covsample <- subset(TEcalling_rates_inner_melt_numbers, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#ratio coverage; re-shape data frame 
#TEcalling_rates_ratio <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_ratio.txt")
#TEcalling_ratio_covsample <- subset(TEcalling_rates_ratio, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#TEcalling_rates_ratio_melt_rates <- melt(TEcalling_rates_ratio, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
#TEcalling_rates_ratio_melt_rates_covsample <- subset(TEcalling_rates_ratio_melt_rates, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#TEcalling_rates_ratio_melt_numbers <- melt(TEcalling_rates_ratio, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_present_calls", "true_absent_calls", "false_pres_calls", "false_absent_calls"))
#TEcalling_rates_ratio_melt_numbers_covsample <- subset(TEcalling_rates_ratio_melt_numbers, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#subsample calling rates 
TEcalling_rates_inner_sub0.20 <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.20.txt")
TEcalling_rates_inner_sub0.20$mean_cov <- 6.5
TEcalling_rates_inner_sub0.32 <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.32.txt")
TEcalling_rates_inner_sub0.32$mean_cov <- 10
TEcalling_rates_inner_sub0.45 <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.45.txt")
TEcalling_rates_inner_sub0.45$mean_cov <- 14
TEcalling_rates_inner_sub0.66 <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.66.txt")
TEcalling_rates_inner_sub0.66$mean_cov <- 20
TEcalling_rates_inner_sub0.83 <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.83.txt")
TEcalling_rates_inner_sub0.83$mean_cov <- 26

#combine subsampled files for graphing 
TEcalling_rates_inner_sub_all <- do.call("rbind", list(TEcalling_rates_inner_sub0.20, TEcalling_rates_inner_sub0.32, TEcalling_rates_inner_sub0.45, TEcalling_rates_inner_sub0.66, TEcalling_rates_inner_sub0.83))
TEcalling_rates_inner_sub_all_covsample <- subset(TEcalling_rates_inner_sub_all, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_inner_sub_all_melt <- melt(TEcalling_rates_inner_sub_all, id.vars = c("population", "cov_limit", "bp_length", "mean_cov"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_inner_sub_all_melt_covsample <- subset(TEcalling_rates_inner_sub_all_melt, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_inner_sub_all_melt_numbers <- melt(TEcalling_rates_inner_sub_all, id.vars = c("population", "cov_limit", "bp_length", "mean_cov"), measure.vars = c("true_present_calls", "true_absent_calls", "false_pres_calls", "false_absent_calls"))
TEcalling_rates_inner_sub_all_melt_numbers_covsample <- subset(TEcalling_rates_inner_sub_all_melt_numbers, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#combine W22 1bp lim calling file
TEcalling_rates_W22_1bplim <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_W22_10bp_1bp.lim.accurate.txt")
TEcalling_rates_W22_1bplim_covsample <- subset(TEcalling_rates_W22_1bplim, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_W22_1bplim_melt <- melt(TEcalling_rates_W22_1bplim, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_W22_1bplim_melt_covsample <- subset(TEcalling_rates_W22_1bplim_melt, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_W22_1bplim_melt_numbers <- melt(TEcalling_rates_W22_1bplim, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_present_calls", "true_absent_calls", "false_pres_calls", "false_absent_calls"))
TEcalling_rates_W22_1bplim_melt_numbers_covsample <- subset(TEcalling_rates_W22_1bplim_melt_numbers, bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#compare W22 normal length and clipped reads 
W22_10bp_calling_inner_normal.clipped <- read_tsv("~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_refB73v4_TEcalling_rates_10bp_normal.clipped.txt")
W22_10bp_calling_inner_normal.clipped_covsample <- subset(W22_10bp_calling_inner_normal.clipped, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
W22_10bp_calling_inner_normal.clipped_melt <- melt(W22_10bp_calling_inner_normal.clipped, id.vars = c("read_len", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
W22_10bp_calling_inner_normal.clipped_melt_covsample <- subset(W22_10bp_calling_inner_normal.clipped_melt, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#read in summary files from 2 category, 4 sides calls 
TEcalling_rates_2cat_4points <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_4points.txt")
TEcalling_rates_2cat_4points_covsample <- subset(TEcalling_rates_2cat_4points,  bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_2cat_4points_melt <- melt(TEcalling_rates_2cat_4points, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_2cat_4points_melt_subset <- subset(TEcalling_rates_2cat_4points_melt, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
#read in summary files from 2 category, 2 side calls 
TEcalling_rates_2cat_inner <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_inner.txt")
TEcalling_rates_2cat_inner_covsample <- subset(TEcalling_rates_2cat_inner,  bp_length != 3 & cov_limit == 1 | bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)
TEcalling_rates_2cat_inner_melt <- melt(TEcalling_rates_2cat_inner, id.vars = c("population", "cov_limit", "bp_length"), measure.vars = c("true_pres_rate_noamb", "true_abs_rate_noamb", "false_pres_rate_noamb", "false_abs_rate_noamb"))
TEcalling_rates_2cat_inner_melt_subset <- subset(TEcalling_rates_2cat_inner_melt, bp_length != 3 & cov_limit == 2 | bp_length != 3 & cov_limit == 5 | bp_length != 3 & cov_limit == 8)

############################################################################
#graphs for presentations - 3 category calls, 4 points and just inner calls 
#n ambiguous calls 
ggplot(TEcalling_rates_4points_covsample, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length") + ylim(40000,135000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_4sides.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_4sides.svg", device = "svg", height = 4, width = 8, units = "in")
ggplot(TEcalling_rates_inner_covsample, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length")+ ylim(40000,135000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_inner.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_inner.svg", device = "svg", height = 4, width = 8, units = "in")
#te calling rates - 4 points and inner coverage 
ggplot(TEcalling_rates_4points_melt_rates_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_covsample.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
ggplot(TEcalling_rates_inner_melt_rates_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype")+ theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_inner_covsample.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_inner_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
#TE calling numbers 
#ggplot(TEcalling_rates_4points_melt_numbers_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Number of calls") + labs(color = "TE fragment length", shape = "Genotype")+ theme(axis.text.x = element_text(angle = 45))
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_4sides_covsample.png", device = "png", height = 4, width = 8, units = "in")
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_4sides_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
#ggplot(TEcalling_rates_inner_melt_numbers_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype")+ theme(axis.text.x = element_text(angle = 45))
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_inner_covsample.png", device = "png", height = 4, width = 8, units = "in")
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_inner_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
#TE calling rates - ratio
#ggplot(TEcalling_rates_ratio_melt_rates_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncratio of coverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_covsample.png", device = "png", height = 4, width = 8, units = "in")
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
#TE calling rates and numbers subsampled files 
ggplot(TEcalling_rates_inner_sub_all_melt_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_subsampled_bams_covsample.png", device = "png", height = 10, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_subsampled_bams_covsample.svg", device = "svg", height = 10, width = 10, units = "in")
#ggplot(TEcalling_rates_inner_sub_all_melt_numbers_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("Number of TEs called in each category\n coverage in inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
#n-ambiguous, subsampled files 
ggplot(TEcalling_rates_inner_sub_all_covsample, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length") 
+ ylim(40000,120000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_subsampled_bams_covsample.png", device = "png", height = 10, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_subsampled_bams_covsample.svg", device = "svg", height = 10, width = 10, units = "in")
#W22 1 bp limit
ggplot(TEcalling_rates_W22_1bplim_covsample, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' in W22\ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length") + ylim(40000,135000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_W22_1bp.lim_covsample.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_W22_1bp.lim_covsample.svg", device = "svg", height = 4, width = 8, units = "in")
ggplot(TEcalling_rates_W22_1bplim_melt_covsample, aes(x = variable, y = value, color = as.factor(bp_length))) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_W22_1bp.lim_covsample.png", device = "png", height = 4, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_W22_1bp.lim_covsample.svg", device = "svg", height = 4, width = 8, units = "in")

#TE calling in W22 normal vs clipped reads 
ggplot(W22_10bp_calling_inner_normal.clipped_melt_covsample, aes(x = variable, y = value, color = as.factor(cov_limit))) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ read_len) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "Coverage limit", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_W22_clipped_v_normal_covsample.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_W22_clipped_v_normal_covsample.svg", device = "svg", height = 6, width = 8, units = "in")
ggplot(W22_10bp_calling_inner_normal.clipped_covsample, aes(x = as.factor(cov_limit), y = n_ambiguous)) + geom_jitter(size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ read_len) + xlab("Coverage Call Threshold") + ylab("Number of Ambiguous Calls") + labs(color = "Coverage limit", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45)) + ylim(40000,120000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_W22_clipped_v_normal_covsample.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_W22_clipped_v_normal_covsample.svg", device = "svg", height = 6, width = 8, units = "in")

#2 category 4 side graphs
ggplot(TEcalling_rates_2cat_4points_melt_subset, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates - only Present/Absent calls \ncoverage on 4 sides of TEs used to call Present/Absent") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_2cat_4sides_covsample.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_2cat_4sides_covsample.svg", device = "svg", height = 6, width = 8, units = "in")
#2 category 2 side graphs
ggplot(TEcalling_rates_2cat_inner_melt_subset, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates - only Present/Absent calls \ncoverage on 2 inner sides of TEs used to call Present/Absent") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_2cat_inner_covsample.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_2cat_inner_covsample.svg", device = "svg", height = 6, width = 8, units = "in")

########################################################################
########################################################################
########################################################################
#Transpose tables for exporting 
#3 category, 4 sides
TEcalling_rates_4points_covsample_t <- t(TEcalling_rates_4points_covsample)
write.table(TEcalling_rates_4points_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_3cat_4points_transpose.txt", sep = "\t")
#3 category, only inner sides of TEs 
TEcalling_rates_inner_covsample_t <- t(TEcalling_rates_inner_covsample)
write.table(TEcalling_rates_inner_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_3cat_inner_transpose.txt", sep = "\t")
#2 category, 4 sides 
TEcalling_rates_2cat_4points_covsample_t <- t(TEcalling_rates_2cat_4points_covsample)
write.table(TEcalling_rates_2cat_4points_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_4points_transpose.txt", sep = "\t")
#2 category, only inner sides of TEs 
TEcalling_rates_2cat_inner_covsample_t <- t(TEcalling_rates_2cat_inner_covsample)
write.table(TEcalling_rates_2cat_inner_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_2cat_inner_transpose.txt", sep = "\t")

#subsample files 
TEcalling_rates_W22_1bplim_covsample_t <- t(TEcalling_rates_W22_1bplim_covsample)
write.table(TEcalling_rates_W22_1bplim_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_W22_1bplim_covsample_transpose.txt", sep = "\t")
TEcalling_rates_inner_sub0.20_t <- t(TEcalling_rates_inner_sub0.20)
write.table(TEcalling_rates_inner_sub0.20_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.20_transpose.txt", sep = "\t")
TEcalling_rates_inner_sub0.32_t <- t(TEcalling_rates_inner_sub0.32)
write.table(TEcalling_rates_inner_sub0.32_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.32_transpose.txt", sep = "\t")
TEcalling_rates_inner_sub0.45_t <- t(TEcalling_rates_inner_sub0.45)
write.table(TEcalling_rates_inner_sub0.45_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/TEcalling_rates_inner_sub0.45_transpose.txt", sep = "\t")

#W22 clipped and normal 
W22_10bp_calling_inner_normal.clipped_covsample_t <- t(W22_10bp_calling_inner_normal.clipped_covsample)
write.table(W22_10bp_calling_inner_normal.clipped_covsample_t, file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/W22_10bp_calling_inner_normal.clipped_covsample_transpose.txt", sep = "\t")

#######################################################################
#######################################################################
#######################################################################
#graphs for supplemental material
ggplot(TEcalling_rates_4points, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_4sides_allsamples.png", device = "png", height = 5, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_4sides_allsamples.svg", device = "svg", height = 5, width = 10, units = "in")
ggplot(TEcalling_rates_inner, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_inner_allsamples.png", device = "png", height = 5, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rate_nambig_inner_allsamples.svg", device = "svg", height = 5, width = 10, units = "in")

#te calling rates 
ggplot(TEcalling_rates_4points_melt_rates, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype")+ theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_allsamples.png", device = "png", height = 5, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_4sides_allsamples.svg", device = "svg", height = 5, width = 10, units = "in")
ggplot(TEcalling_rates_inner_melt_rates, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_inner_allsamples.png", device = "png", height = 5, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_inner_allsamples.svg", device = "svg", height = 5, width = 10, units = "in")

#subsample files 

#TE calling rates and numbers subsampled files 
ggplot(TEcalling_rates_inner_sub_all_melt, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_subsampled_bams.png", device = "png", height = 10, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_rates_innner_subsampled_bams.svg", device = "svg", height = 10, width = 10, units = "in")
#ggplot(TEcalling_rates_inner_sub_all_melt_numbers_covsample, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("Number of TEs called in each category\n coverage in inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
#n-ambiguous, subsampled files 
ggplot(TEcalling_rates_inner_sub_all, aes(x = population, y = n_ambiguous, color = as.factor(bp_length))) + geom_point(size = 3) + ggtitle("Number of TEs called as 'Ambiguous' \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(mean_cov ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + ylab("Number Ambiguous TE calls") + xlab("Genotype") + labs(color = "TE fragment length") + ylim(40000,120000)
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_subsampled_bams.png", device = "png", height = 10, width = 10, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_nambig_innner_subsampled_bams.svg", device = "svg", height = 10, width = 10, units = "in")
#TE calling in W22 normal vs clipped reads 
ggplot(W22_10bp_calling_inner_normal.clipped_melt, aes(x = variable, y = value, color = as.factor(cov_limit))) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ read_len) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + xlab("Call Type") + ylab("Call rate") + labs(color = "Coverage limit", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_W22_clipped_v_normal.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_W22_clipped_v_normal.svg", device = "svg", height = 6, width = 8, units = "in")

#TE calling numbers 
#ggplot(TEcalling_rates_4points_melt_numbers, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside and outside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + xlab("Call Type") + ylab("Number of calls") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_4sides_allsamples.png", device = "png", height = 4, width = 8, units = "in")
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_4sides_allsamples.svg", device = "svg", height = 4, width = 8, units = "in")
#ggplot(TEcalling_rates_inner_melt_numbers, aes(x = variable, y = value, color = as.factor(bp_length), shape = population)) + geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.75) + ggtitle("TE accuracy calling rates \ncoverage on inside of TEs used to call Present/Absent/Ambiguous") + facet_grid(. ~ cov_limit) + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) + xlab("Call Type") + ylab("Call rate") + labs(color = "TE fragment length", shape = "Genotype") + theme(axis.text.x = element_text(angle = 45))
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_inner_allsamples.png", device = "png", height = 4, width = 8, units = "in")
#ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/TEcalling_numbers_inner_allsamples.svg", device = "svg", height = 4, width = 8, units = "in")

#alignment coverage distribution
group2_B73_cov <- read_tsv(file = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/B73_group2_avg_cov_try2.txt")
ggplot(group2_B73_cov, aes(x = coverage)) + geom_histogram(bins = 40) + ggtitle("Genome wide mean coverage of 430 WiDiv genotypes") + xlab("Mean coverage")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/B73_group2_avg_cov_distribution.png", device = "png", height = 6, width = 8, units = "in")
ggsave(filename = "~/Documents/TE_project/true_false_pos_calling/bedtools_true_false_pos/figures/B73_group2_avg_cov_distribution.svg", device = "svg", height = 6, width = 8, units = "in")

#end of script 


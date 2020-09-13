#calculate calling accuracy of bedtools multicov for calling TEs as present or absent 
library(tidyverse)

#read in files 
Mo17_refB73v4_TEcov_1bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_TEcov_1bp_Mo17status.txt")
Mo17_refB73v4_TEcov_3bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_TEcov_3bp_Mo17status.txt")
Mo17_refB73v4_TEcov_5bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_TEcov_5bp_Mo17status.txt")
Mo17_refB73v4_TEcov_10bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_TEcov_10bp_Mo17status.txt")
PH207_refB73v4_TEcov_1bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73v4_TEcov_1bp_PH207status.txt")
PH207_refB73v4_TEcov_3bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73v4_TEcov_3bp_PH207status.txt")
PH207_refB73v4_TEcov_5bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73v4_TEcov_5bp_PH207status.txt")
PH207_refB73v4_TEcov_10bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73v4_TEcov_10bp_PH207status.txt")
W22_refB73v4_TEcov_1bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73v4_TEcov_1bp_W22status.txt")
W22_refB73v4_TEcov_3bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73v4_TEcov_3bp_W22status.txt")
W22_refB73v4_TEcov_5bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73v4_TEcov_5bp_W22status.txt")
W22_refB73v4_TEcov_10bp <- read_tsv(file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73v4_TEcov_10bp_W22status.txt")

#functions to call calls as true or false positives 
calc_sensitivity_cov_counts <- function(dataframe_count) {
        dataframe_count$calls_2bp_lim <- as.character("NA")
        dataframe_count$calls_3bp_lim <- as.character("NA")
        dataframe_count$calls_4bp_lim <- as.character("NA")
        dataframe_count$calls_5bp_lim <- as.character("NA")
        dataframe_count$calls_6bp_lim <- as.character("NA")
        dataframe_count$calls_7bp_lim <- as.character("NA")
        dataframe_count$calls_8bp_lim <- as.character("NA")
        cov_limits <- c(2,3,4,5,6,7,8)
        for (num in cov_limits) {
        	abs = num - 1
        	for (i in 1:nrow(dataframe_count)) {
                if (dataframe_count[i,2] == 'Present' & dataframe_count[i,3] >= num & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num & dataframe_count[i,6] >= num) {
                        dataframe_count[i,5+num] <- "call_present" #true positive
                } else if (dataframe_count[i,2] == 'Absent' & dataframe_count[i,3] < abs & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs & dataframe_count[i,6] < abs) {
                        dataframe_count[i,5+num] <- "call_absent" #true negative
                } else if (dataframe_count[i,2] == 'Absent.SD' & dataframe_count[i,3] < abs & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs & dataframe_count[i,6] < abs) {
                        dataframe_count[i,5+num] <- "call_absent" #true negative
                } else if (dataframe_count[i,2] == 'Present' & dataframe_count[i,3] < abs & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs & dataframe_count[i,6] < abs) {
                        dataframe_count[i,5+num] <- "call_false_absent" #false negative, should be present but coverage is too low
                } else if (dataframe_count[i,2] == 'Absent' & dataframe_count[i,3] >= num & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num & dataframe_count[i,6] >= num) {
                        dataframe_count[i,5+num] <- "call_false_present" #false positive, should be absent but has high coverage
                } else if (dataframe_count[i,2] == 'Absent.SD' & dataframe_count[i,3] >= num & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num & dataframe_count[i,6] >= num) {
                        dataframe_count[i,5+num] <- "call_false_present" #false positive, should be absent but has high coverage
                } else if (dataframe_count[i,2] == 'Unresolved') {
                		dataframe_count[i,5+num] <- "undefined"
                }  else {
                		dataframe_count[i,5+num] <- "ambiguous"
                	}
            }
	}
        return(dataframe_count)
}

calc_sensitivity_cov_counts_innercov <- function(dataframe_count) {
	
        dataframe_count$calls_2bp_lim <- as.character("NA")
        dataframe_count$calls_3bp_lim <- as.character("NA")
        dataframe_count$calls_4bp_lim <- as.character("NA")
        dataframe_count$calls_5bp_lim <- as.character("NA")
        dataframe_count$calls_6bp_lim <- as.character("NA")
        dataframe_count$calls_7bp_lim <- as.character("NA")
        dataframe_count$calls_8bp_lim <- as.character("NA")
        cov_limits <- c(2,3,4,5,6,7,8)
        for (num in cov_limits) {
        	abs = num - 1
        	for (i in 1:nrow(dataframe_count)) {
                if (dataframe_count[i,2] == 'Present' & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num) {
                        dataframe_count[i,5+num] <- "call_present" #true positive
                } else if (dataframe_count[i,2] == 'Absent' & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs) {
                        dataframe_count[i,5+num] <- "call_absent" #true negative
                } else if (dataframe_count[i,2] == 'Absent.SD' & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs) {
                        dataframe_count[i,5+num] <- "call_absent" #true negative
                } else if (dataframe_count[i,2] == 'Present' & dataframe_count[i,4] < abs & dataframe_count[i,5] < abs) {
                        dataframe_count[i,5+num] <- "call_false_absent" #false negative, should be present but coverage is too low
                } else if (dataframe_count[i,2] == 'Absent' & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num) {
                        dataframe_count[i,5+num] <- "call_false_present" #false positive, should be absent but has high coverage
                } else if (dataframe_count[i,2] == 'Absent.SD' & dataframe_count[i,4] >= num & dataframe_count[i,5] >= num) {
                        dataframe_count[i,5+num] <- "call_false_present" #false positive, should be absent but has high coverage
                } else if (dataframe_count[i,2] == 'Unresolved') {
                		dataframe_count[i,5+num] <- "undefined"
                }  else {
                		dataframe_count[i,5+num] <- "ambiguous"
                	}
            }
	}
        return(dataframe_count)
}


Mo17_refB73v4_TEcov_1bp_calls <- calc_sensitivity_cov_counts(Mo17_refB73v4_TEcov_1bp)
write.table(Mo17_refB73v4_TEcov_1bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_1bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_3bp_calls <- calc_sensitivity_cov_counts(Mo17_refB73v4_TEcov_3bp)
write.table(Mo17_refB73v4_TEcov_3bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_3bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_5bp_calls <- calc_sensitivity_cov_counts(Mo17_refB73v4_TEcov_5bp)
write.table(Mo17_refB73v4_TEcov_5bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_5bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_10bp_calls <- calc_sensitivity_cov_counts(Mo17_refB73v4_TEcov_10bp)
write.table(Mo17_refB73v4_TEcov_10bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_10bp.tf.pos.4sides.txt", sep = "\t", row.names = F)

PH207_refB73v4_TEcov_1bp_calls <- calc_sensitivity_cov_counts(PH207_refB73v4_TEcov_1bp)
write.table(PH207_refB73v4_TEcov_1bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_1bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_3bp_calls <- calc_sensitivity_cov_counts(PH207_refB73v4_TEcov_3bp)
write.table(PH207_refB73v4_TEcov_3bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_3bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_5bp_calls <- calc_sensitivity_cov_counts(PH207_refB73v4_TEcov_5bp)
write.table(PH207_refB73v4_TEcov_5bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_5bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_10bp_calls <- calc_sensitivity_cov_counts(PH207_refB73v4_TEcov_10bp)
write.table(PH207_refB73v4_TEcov_10bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_10bp.tf.pos.4sides.txt", sep = "\t", row.names = F)

W22_refB73v4_TEcov_1bp_calls <- calc_sensitivity_cov_counts(W22_refB73v4_TEcov_1bp)
write.table(W22_refB73v4_TEcov_1bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_1bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_3bp_calls <- calc_sensitivity_cov_counts(W22_refB73v4_TEcov_3bp)
write.table(W22_refB73v4_TEcov_3bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_3bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_5bp_calls <- calc_sensitivity_cov_counts(W22_refB73v4_TEcov_5bp)
write.table(W22_refB73v4_TEcov_5bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_5bp.tf.pos.4sides.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_10bp_calls <- calc_sensitivity_cov_counts(W22_refB73v4_TEcov_10bp)
write.table(W22_refB73v4_TEcov_10bp_calls, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_10bp.tf.pos.4sides.txt", sep = "\t", row.names = F)

Mo17_refB73v4_TEcov_1bp_calls_inner <- calc_sensitivity_cov_counts_innercov(Mo17_refB73v4_TEcov_1bp)
write.table(Mo17_refB73v4_TEcov_1bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_1bp.tf.pos_inner.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_3bp_calls_inner <- calc_sensitivity_cov_counts_innercov(Mo17_refB73v4_TEcov_3bp)
write.table(Mo17_refB73v4_TEcov_3bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_3bp.tf.pos_inner.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_5bp_calls_inner <- calc_sensitivity_cov_counts_innercov(Mo17_refB73v4_TEcov_5bp)
write.table(Mo17_refB73v4_TEcov_5bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_5bp.tf.pos_inner.txt", sep = "\t", row.names = F)
Mo17_refB73v4_TEcov_10bp_calls_inner <- calc_sensitivity_cov_counts_innercov(Mo17_refB73v4_TEcov_10bp)
write.table(Mo17_refB73v4_TEcov_10bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73_TEcov_10bp.tf.pos_inner.txt", sep = "\t", row.names = F)

PH207_refB73v4_TEcov_1bp_calls_inner <- calc_sensitivity_cov_counts_innercov(PH207_refB73v4_TEcov_1bp)
write.table(PH207_refB73v4_TEcov_1bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_1bp.tf.pos_inner.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_3bp_calls_inner <- calc_sensitivity_cov_counts_innercov(PH207_refB73v4_TEcov_3bp)
write.table(PH207_refB73v4_TEcov_3bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_3bp.tf.pos_inner.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_5bp_calls_inner <- calc_sensitivity_cov_counts_innercov(PH207_refB73v4_TEcov_5bp)
write.table(PH207_refB73v4_TEcov_5bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_5bp.tf.pos_inner.txt", sep = "\t", row.names = F)
PH207_refB73v4_TEcov_10bp_calls_inner <- calc_sensitivity_cov_counts_innercov(PH207_refB73v4_TEcov_10bp)
write.table(PH207_refB73v4_TEcov_10bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/PH207_refB73_TEcov_10bp.tf.pos_inner.txt", sep = "\t", row.names = F)

W22_refB73v4_TEcov_1bp_calls_inner <- calc_sensitivity_cov_counts_innercov(W22_refB73v4_TEcov_1bp)
write.table(W22_refB73v4_TEcov_1bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_1bp.tf.pos_inner.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_3bp_calls_inner <- calc_sensitivity_cov_counts_innercov(W22_refB73v4_TEcov_3bp)
write.table(W22_refB73v4_TEcov_3bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_3bp.tf.pos_inner.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_5bp_calls_inner <- calc_sensitivity_cov_counts_innercov(W22_refB73v4_TEcov_5bp)
write.table(W22_refB73v4_TEcov_5bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_5bp.tf.pos_inner.txt", sep = "\t", row.names = F)
W22_refB73v4_TEcov_10bp_calls_inner <- calc_sensitivity_cov_counts_innercov(W22_refB73v4_TEcov_10bp)
write.table(W22_refB73v4_TEcov_10bp_calls_inner, file = "/home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/W22_refB73_TEcov_10bp.tf.pos_inner.txt", sep = "\t", row.names = F)

#end of script 


# joining matrix
# to join matrix into a non-redundant matrix,
# reformt the bed file from the TE contrast file 

# B73 ref 
B73_bed <- read.csv("B73_TE_fmt.bed", sep = "\t", header= FALSE)
colnames(B73_bed) <- c("TE_ID","B73.coordinates")
refB73_freq <- read.csv("refB73_widiv_freq.csv", header= TRUE)
head(refB73_freq)
refB73_matrix <- left_join(refB73_freq[2:6],B73_bed)

# Mo17 ref
Mo17_bed <- read.csv("Mo17_TE_fmt.bed", sep = "\t", header= FALSE)
colnames(Mo17_bed) <- c("TE_ID","Mo17.coordinates")
refMo17_freq <- read.csv("refMo17_widiv_freq.csv", header= TRUE)
head(refMo17_freq)
refMo17_matrix <- left_join(refMo17_freq[2:6],Mo17_bed)

# PH207 ref 
PH207_bed <- read.csv("PH207_TE_fmt.bed", sep = "\t", header= FALSE)
colnames(PH207_bed) <- c("TE_ID","PH207.coordinates")
refPH207_freq <- read.csv("refPH207_widiv_freq.csv", header= TRUE)
head(refPH207_freq)
refPH207_matrix <- left_join(refPH207_freq[2:6],PH207_bed)

# W22
W22_bed <- read.csv("W22_TE_fmt.bed", sep = "\t", header= FALSE)
colnames(W22_bed) <- c("TE_ID","W22.coordinates")
refW22_freq <- read.csv("refW22_widiv_freq.csv", header= TRUE)
head(refW22_freq)
refW22_matrix <- left_join(refW22_freq[2:6],W22_bed)

sum(is.na(refW22_matrix$W22.coordinates))
sum(is.na(refPH207_matrix$PH207.coordinates))
sum(is.na(refMo17_matrix$Mo17.coordinates))
sum(is.na(refB73_matrix$B73.coordinates))

head(refB73_matrix)

# base for TE frequency join
base <- read.csv("base_for_non_TE_joint.txt", sep="\t",header=TRUE)
# prepare matrix for join

refB73_matrix_for_join <- refB73_matrix %>% filter(B73.coordinates != "NA") %>% 
  select("B73.coordinates","TE_ID","pop_freq")
refMo17_matrix_for_join <- refMo17_matrix %>% filter(Mo17.coordinates != "NA") %>% 
  select("Mo17.coordinates","TE_ID","pop_freq")
refPH207_matrix_for_join <- refPH207_matrix %>% filter(PH207.coordinates != "NA") %>% 
  select("PH207.coordinates","TE_ID","pop_freq")
refW22_matrix_for_join <- refW22_matrix %>% filter(W22.coordinates != "NA") %>% 
  select("W22.coordinates","TE_ID","pop_freq")


freq_matrix <- left_join(base,refB73_matrix_for_join,by="B73.coordinates") %>% left_join(refMo17_matrix_for_join,by="Mo17.coordinates") %>% 
  left_join(refPH207_matrix_for_join, by="PH207.coordinates") %>% left_join(refW22_matrix_for_join, by="W22.coordinates")
colnames(freq_matrix) <- c("TE","disjoined","B73.coordinates","W22.coordinates","Mo17.coordinates","PH207.coordinates","source","present.genos","family", "variable","superfam","B73_TE_ID","pop_freq_B73","TE_ID_Mo17","pop_freq_Mo17","TE_ID_PH207","pop_freq_PH207","TE_ID_W22","pop_freq_W22")
freq_matrix$prop_present <- rowMeans(subset(freq_matrix, select = c(pop_freq_B73, pop_freq_Mo17,pop_freq_PH207,pop_freq_W22)), na.rm = TRUE)

# LTR similarity and nest info 
LTR_age <- read.csv("allTEs_LTRsimilarity.txt", sep="\t",header = FALSE)
colnames(LTR_age) <- c("TE","LTR_age")
# add in nest status info
nest_stats <- read.csv("TE_nest_info.txt",sep = '\t', header =FALSE)
colnames(nest_stats) <- c("TE","nested_status")

matrix_for_analysis = freq_matrix %>% filter(prop_present != "NaN") %>% 
  select("TE","disjoined","source","family", "variable","superfam","prop_present") %>% 
  left_join(LTR_age,by="TE") %>% left_join(nest_stats,by="TE")

# add in TE_order
matrix_for_analysis$order <-ifelse(matrix_for_analysis$superfam == "DHH",rr2 <-"Helitron",
                                   ifelse(matrix_for_analysis$superfam == "DTA",rr2 <-"TIR",
                                          ifelse(matrix_for_analysis$superfam == "DTC", rr2<-"TIR",
                                                 ifelse(matrix_for_analysis$superfam == "DTH", rr2<-"TIR",
                                                        ifelse(matrix_for_analysis$superfam == "DTM", rr2<-"TIR",
                                                               ifelse(matrix_for_analysis$superfam == "DTT", rr2<-"TIR",
                                                                      ifelse(matrix_for_analysis$superfam == "DTX", rr2<-"TIR",
                                                                             ifelse(matrix_for_analysis$superfam == "RIL", rr2<-"LINE",
                                                                                    ifelse(matrix_for_analysis$superfam == "RIT", rr2<-"LINE",
                                                                                           ifelse(matrix_for_analysis$superfam == "RLC", rr2<-"LTR",
                                                                                                  ifelse(matrix_for_analysis$superfam == "RLG", rr2<-"LTR",
                                                                                                         ifelse(matrix_for_analysis$superfam == "RLX", rr2<-"LTR",
                                                                                                                rr2<-"SINE"))))))))))))

# final meta info 
meta_info <- read.csv("all_TE_meta.txt",sep = '\t',header=TRUE)
colnames(meta_info) <- c("TE","TE_len","class","fam_size","genomic_loc")
final_matrix <- left_join(matrix_for_analysis,meta_info) 
write.csv(final_matrix,"Widiv_TE_variation_matrix_YQ.csv")

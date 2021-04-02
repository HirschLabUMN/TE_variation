# Only TEs that resequencing mapping to the ref itself can call out were kept for the analysis 

library(tidyverse)
# refB73 

refB73_conver_PAV = read.csv("refB73v4_order_widiv_509.csv_present_probability_wide.csv",header=TRUE)
refB73_conver_PAV[refB73_conver_PAV >= 0.7 ] <- "present"
refB73_conver_PAV[refB73_conver_PAV <= 0.3 ] <- "absent"

refB73_conver_PAV$present_count <- rowSums(refB73_conver_PAV[2:510] == "present")
refB73_conver_PAV$absent_count <- rowSums(refB73_conver_PAV[2:510] == "absent")
refB73_conver_PAV$ambig <- 509- refB73_conver_PAV$present_count - refB73_conver_PAV$absent_count
refB73_conver_PAV$pop_freq = refB73_conver_PAV$present_count/(refB73_conver_PAV$present_count+refB73_conver_PAV$absent_count)
refB73_conver_PAV_filter = refB73_conver_PAV %>% filter(B73 == "present") 
refB73_conver_PAV_filter %>% select("TE_ID","present_count","absent_count","pop_freq","ambig") %>% write.csv("refB73_widiv_freq.csv")

# refMo17 

refMo17_conver_PAV = read.csv("refMo17_order_widiv_509.csv_present_probability_wide.csv",header=TRUE)
refMo17_conver_PAV[refMo17_conver_PAV >= 0.7 ] <- "present"
refMo17_conver_PAV[refMo17_conver_PAV <= 0.3 ] <- "absent"

refMo17_conver_PAV$present_count <- rowSums(refMo17_conver_PAV[2:510] == "present")
refMo17_conver_PAV$absent_count <- rowSums(refMo17_conver_PAV[2:510] == "absent")
refMo17_conver_PAV$ambig <- 509- refMo17_conver_PAV$present_count - refMo17_conver_PAV$absent_count
refMo17_conver_PAV$pop_freq = refMo17_conver_PAV$present_count/(refMo17_conver_PAV$present_count+refMo17_conver_PAV$absent_count)
refMo17_conver_PAV_filter = refMo17_conver_PAV %>% filter(Mo17 == "present") 
refMo17_conver_PAV_filter %>% select("TE_ID","present_count","absent_count","pop_freq","ambig") %>% write.csv("refMo17_widiv_freq.csv")

# refPH207 

refPH207_conver_PAV = read.csv("refPH207_order_widiv_509.csv_present_probability_wide.csv",header=TRUE)
refPH207_conver_PAV[refPH207_conver_PAV >= 0.7 ] <- "present"
refPH207_conver_PAV[refPH207_conver_PAV <= 0.3 ] <- "absent"

refPH207_conver_PAV$present_count <- rowSums(refPH207_conver_PAV[2:510] == "present")
refPH207_conver_PAV$absent_count <- rowSums(refPH207_conver_PAV[2:510] == "absent")
refPH207_conver_PAV$ambig <- 509- refPH207_conver_PAV$present_count - refPH207_conver_PAV$absent_count
refPH207_conver_PAV$pop_freq = refPH207_conver_PAV$present_count/(refPH207_conver_PAV$present_count+refPH207_conver_PAV$absent_count)
refPH207_conver_PAV_filter = refPH207_conver_PAV %>% filter(PH207 == "present") 
refPH207_conver_PAV_filter %>% select("TE_ID","present_count","absent_count","pop_freq","ambig") %>% write.csv("refPH207_widiv_freq.csv")

# refW22 

refW22_conver_PAV = read.csv("refW22_order_widiv_509.csv_present_probability_wide.csv",header=TRUE)
refW22_conver_PAV[refW22_conver_PAV >= 0.7 ] <- "present"
refW22_conver_PAV[refW22_conver_PAV <= 0.3 ] <- "absent"

refW22_conver_PAV$present_count <- rowSums(refW22_conver_PAV[2:510] == "present")
refW22_conver_PAV$absent_count <- rowSums(refW22_conver_PAV[2:510] == "absent")
refW22_conver_PAV$ambig <- 509- refW22_conver_PAV$present_count - refW22_conver_PAV$absent_count
refW22_conver_PAV$pop_freq = refW22_conver_PAV$present_count/(refW22_conver_PAV$present_count+refW22_conver_PAV$absent_count)
refW22_conver_PAV_filter = refW22_conver_PAV %>% filter(W22 == "present") 
refW22_conver_PAV_filter %>% select("TE_ID","present_count","absent_count","pop_freq","ambig") %>% write.csv("refW22_widiv_freq.csv")

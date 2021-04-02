library(data.table)
library(tidyverse)
files <- list.files(path="/home/hirschc1/qiuxx221/TE_variation/rf_prediction/prediction/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  df <- fread(x, header = FALSE, data.table = FALSE)
  colnames(df) <- c("TE_ID","start","end","Sequenced_genome","TE_fam","TE_order","absent","present")
  df_select = df %>% select(TE_ID,Sequenced_genome,present) 
  df_select %>% 
    spread(Sequenced_genome, present) %>% write.csv(paste0(x, "_present_probability_wide.csv"), row.names = FALSE)
})

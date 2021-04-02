setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
TE_matrix <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv",header = TRUE)

TE_matrix %>% filter(LTR_age != "NA") %>% 
  ggplot(aes(x=LTR_age, color=nested_status)) + 
  geom_density() +  theme_classic() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),labels = c("Nested", "Non-nested")) +
  xlab("LTR Similarity (%)") + ylab("Density") + labs(name = "TE Nested Status") +
  theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
                                                   axis.text.y = element_text(size=12), 
                                                   axis.text.x = element_text(size=12),
                                                   legend.text = element_text(size=12),
                                                   legend.title=element_blank()) +
  theme(legend.position = c(.2,.85)) + expand_limits(x = 84)+ expand_limits(x = 101)

dim(TE_matrix %>% filter(LTR_age != "NA"))

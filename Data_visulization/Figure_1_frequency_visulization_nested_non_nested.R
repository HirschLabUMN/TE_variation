setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
TE_matrix <- read.csv("Widiv_TE_variation_matrix_revision_v1_fmt.csv",header = TRUE)

library(tidyverse)
# nested TEs 
table(TE_matrix$nested_status)
# nested_TE non_nested_TE 
# 253714        191704 

TE_matrix$class_stack<-ifelse(TE_matrix$prop_present <= 0.2,rr2<-"0-20%",
                                        ifelse(TE_matrix$prop_present>0.2 &TE_matrix$prop_present<=0.8,rr2<-"20%-80%",
                                               ifelse(TE_matrix$prop_present>0.8 & TE_matrix$prop_present<=1, rr2<-"80%-100%",
                                                      rr2<-"")))

TE_matrix$age_stack<-ifelse(TE_matrix$LTR_age <= 95,rr2<-"Old",
                                      ifelse(TE_matrix$LTR_age > 95 &TE_matrix$LTR_age<= 99,rr2<-"Young",
                                             ifelse(TE_matrix$LTR_age>99 & TE_matrix$LTR_age<=100, rr2<-"Very Young",
                                                    rr2<-"")))


# remove LINE and SINE for the visulization 
three_class_TE_matrix = TE_matrix %>% 
  filter(!grepl('SINE', order)) %>% filter(!grepl('LINE', order))

three_class_TE_matrix$order <- factor(three_class_TE_matrix$order,levels = c("LTR","Helitron","TIR"))

## Nested TE group 
nested_TE <- three_class_TE_matrix %>% filter(nested_status == "nested_TE") 

table(nested_TE$order)
#      LTR Helitron      TIR 
# 92538    13589   146761

Nested_TE_panel = ggplot(nested_TE, aes(x = prop_present)) + geom_histogram() +
  geom_vline(data = nested_TE, mapping = aes(xintercept = 0.2),linetype = "dashed",color="blue") + geom_vline(data = nested_TE, mapping = aes(xintercept = 0.8),linetype = "dashed",color="blue") + 
  facet_wrap(nested_TE$order,scales = "free",strip.position = "bottom") + scale_x_continuous(name="Population Frequency") + scale_y_continuous("Number of Nested TEs")+theme_classic() + 
  theme(strip.background =element_rect(fill="grey"),text = element_text(size = 12)) 





# Non-nested TE group
non_nested_TE <- three_class_TE_matrix %>% filter(nested_status == "non_nested_TE") 
table(non_nested_TE$order)
# LTR Helitron      TIR 
# 97835    13791    78560 

Non_nested_TE_panel = ggplot(non_nested_TE, aes(x = prop_present)) + geom_histogram() +
  geom_vline(data = non_nested_TE, mapping = aes(xintercept = 0.2),linetype = "dashed",color="blue") + geom_vline(data = non_nested_TE, mapping = aes(xintercept = 0.8),linetype = "dashed",color="blue") + 
  facet_wrap(non_nested_TE$order,scales = "free",strip.position = "bottom") + scale_x_continuous(name="Population Frequency") + scale_y_continuous("Number of Non-Nested TEs")+theme_classic() + 
  theme(strip.background =element_rect(fill="grey"),text = element_text(size = 12)) 

grid.arrange(Nested_TE_panel,Non_nested_TE_panel,nrow = 2)




# getting stats for each group by every 20% 
options(pillar.sigfigs=100)
nested_TE %>% select(order,prop_present,class_stack) %>% count(order, class_stack, sort = TRUE) %>% 
  group_by(order) %>% mutate(percent = round(n/sum(n),4))

non_nested_TE %>% select(order,prop_present,class_stack) %>% count(order, class_stack, sort = TRUE) %>% 
  group_by(order) %>% mutate(percent = round(100*n/sum(n),4))


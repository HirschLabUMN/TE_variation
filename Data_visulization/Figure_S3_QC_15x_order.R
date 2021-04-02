setwd("~/Desktop/TE paper/Revision/caret_rf/Final/15x/")
library(tidyverse)
library(caret)
library(lemon)
PH207_tester_15x <- read.csv(file = "PH207_tester_15x.csv")
Mo17_tester_15x <- read.csv(file="Mo17_tester_15x.csv")
B73_tester_15x <- read.csv(file="B73_tester_15x.csv")
W22_tester_15x <- read.csv(file="W22_tester_15x.csv")


B73_training_model_15x_order <- readRDS("no_B73_training_15_default_mtry123_order.rds")
W22_training_model_15x_order <- readRDS("no_W22_training_15_default_mtry123_order.rds")
Mo17_training_model_15x_order <- readRDS("no_Mo17_training_15_default_mtry123_order.rds")
PH207_training_model_15x_order <- readRDS("no_PH207_training_15_default_mtry123_order.rds")

predict_prob_B73_15x_order <- predict(B73_training_model_15x_order, newdata=B73_tester_15x, type="prob")
B73_pred_matrix_15x_order <- cbind(B73_tester_15x,predict_prob_B73_15x_order)

B73_seq_on_no_B73_model_present_15x_order <- ggplot(data = B73_pred_matrix_15x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()
B73_seq_on_no_B73_model_absent_15x_order <- ggplot(data = B73_pred_matrix_15x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()

pdf("15x_B73_reseq_on_no_B73_model_15x_order.pdf",width=8, height=5)
grid.arrange(B73_seq_on_no_B73_model_present_15x_order,B73_seq_on_no_B73_model_absent_15x_order,nrow=1)
dev.off()

ref_genome =c("Mo17","PH207","W22v12")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(B73_pred_matrix_15x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = B73_pred_matrix_15x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = B73_pred_matrix_15x_order %>% filter(absent >=j & Ref_genome ==i)
    total_classifed = as.data.frame(dim(pred_present))[1,] + as.data.frame(dim(pred_absent))[1,]
    present_call = as.data.frame(table(pred_present$TE_PAV))
    false_pos = present_call[1,2]/(present_call[1,2]+present_call[2,2])
    true_pos = present_call[2,2]/(present_call[1,2]+present_call[2,2])
    absent_call= as.data.frame(table(pred_absent$TE_PAV))
    false_neg = absent_call[2,2]/(absent_call[1,2]+absent_call[2,2])
    true_neg = absent_call[1,2]/(absent_call[1,2]+absent_call[2,2])
    pct_captured=total_classifed/total_count 
    QC_stats= as.data.frame(rbind(true_pos,false_pos,true_neg,false_neg,pct_captured))
    QC_stats$Prob_cutoff <- j
    QC_stats$Ref_genome <- i
    colnames(QC_stats[1]) <- "rate"
    QC_stats$metric <- c("true_pos","false_pos","true_neg","false_neg","pct_captured") 
    print(QC_stats)
    out=rbind(out,QC_stats)
  }
}

as.data.frame(out) %>% write.csv("QC_B73_reseq_on_no_B73_model_matrix_15x_order.csv")
visualize_QC_B73_no_B73_model_15x_order <- read.csv(file ="QC_B73_reseq_on_no_B73_model_matrix_15x_order.csv", header= TRUE)
# QC_B73_reseq_on_no_B73_model_matrix_15x_order.csv
visualize_QC_B73_no_B73_model_15x_order$metric <- factor(visualize_QC_B73_no_B73_model_15x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))
B73_reseq_QC_15x_order = visualize_QC_B73_no_B73_model_15x_order %>% 
  ggplot(aes(Prob_cutoff,V1,colour = metric))+ facet_wrap(~Ref_genome) + 
  geom_point(position = position_dodge(width = 0.05)) + 
  xlab("Probability Cutoff") +ylab ("Percentage") + 
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#332288","#CC79A7"), 
                     name="Metrics",
                     breaks=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"),
                     labels=c("False Negative", "False Positive","True Negative","True Positive","Total Pecentage of TE Captured")) + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=12),
        legend.text = element_text(size=12),plot.title = element_text(hjust = 0.5)) + ggtitle("B73 Resequencing Reads") + 
  theme(axis.title.x=element_blank())



######### Mo17 ############
# using no_Mo17 training model to predict Mo17 resequencing reads
predict_prob_Mo17_15x_order <- predict(Mo17_training_model_15x_order, newdata=Mo17_tester_15x, type="prob")
Mo17_pred_matrix_15x_order <- cbind(Mo17_tester_15x,predict_prob_Mo17_15x_order)

Mo17_seq_on_no_Mo17_model_present_15x_order <- ggplot(data = Mo17_pred_matrix_15x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()
Mo17_seq_on_no_Mo17_model_absent_15x_order <- ggplot(data = Mo17_pred_matrix_15x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()

pdf("15x_Mo17_reseq_on_no_Mo17_model_15x_order.pdf",width=8, height=5)
grid.arrange(Mo17_seq_on_no_Mo17_model_present_15x_order,Mo17_seq_on_no_Mo17_model_absent_15x_order,nrow=1)
dev.off()

ref_genome =c("B73v4","PH207","W22v12")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(Mo17_pred_matrix_15x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = Mo17_pred_matrix_15x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = Mo17_pred_matrix_15x_order %>% filter(absent >=j & Ref_genome ==i)
    total_classifed = as.data.frame(dim(pred_present))[1,] + as.data.frame(dim(pred_absent))[1,]
    present_call = as.data.frame(table(pred_present$TE_PAV))
    false_pos = present_call[1,2]/(present_call[1,2]+present_call[2,2])
    true_pos = present_call[2,2]/(present_call[1,2]+present_call[2,2])
    absent_call= as.data.frame(table(pred_absent$TE_PAV))
    false_neg = absent_call[2,2]/(absent_call[1,2]+absent_call[2,2])
    true_neg = absent_call[1,2]/(absent_call[1,2]+absent_call[2,2])
    pct_captured=total_classifed/total_count 
    QC_stats= as.data.frame(rbind(true_pos,false_pos,true_neg,false_neg,pct_captured))
    QC_stats$Prob_cutoff <- j
    QC_stats$Ref_genome <- i
    colnames(QC_stats[1]) <- "rate"
    QC_stats$metric <- c("true_pos","false_pos","true_neg","false_neg","pct_captured") 
    print(QC_stats)
    out=rbind(out,QC_stats)
  }
}
as.data.frame(out) %>% write.csv("QC_Mo17_reseq_on_no_Mo17_model_matrix_15x_order.csv")
visualize_QC_Mo17_no_Mo17_model_15x_order <- read.csv(file ="QC_Mo17_reseq_on_no_Mo17_model_matrix_15x_order.csv", header= TRUE)
visualize_QC_Mo17_no_Mo17_model_15x_order$metric <- factor(visualize_QC_Mo17_no_Mo17_model_15x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))
Mo17_reseq_QC_15x_order = visualize_QC_Mo17_no_Mo17_model_15x_order %>% 
  ggplot(aes(Prob_cutoff,V1,colour = metric))+ facet_wrap(~Ref_genome) + 
  geom_point(position = position_dodge(width = 0.05)) + 
  xlab("Probability Cutoff") +ylab ("Percentage") + 
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#332288","#CC79A7"), 
                     name="Metrics",
                     breaks=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"),
                     labels=c("False Negative", "False Positive","True Negative","True Positive","Total Pecentage of TE Captured")) + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=12),
        legend.text = element_text(size=12),plot.title = element_text(hjust = 0.5)) + ggtitle("Mo17 Resequencing Reads") + 
  theme(axis.title.x=element_blank())



######### PH207 ############
# using no_PH207 training model to predict PH207 resequencing reads
predict_prob_PH207_15x_order <- predict(PH207_training_model_15x_order, newdata=PH207_tester_15x, type="prob")
PH207_pred_matrix_15x_order <- cbind(PH207_tester_15x,predict_prob_PH207_15x_order)

PH207_seq_on_no_PH207_model_present_15x_order <- ggplot(data = PH207_pred_matrix_15x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()
PH207_seq_on_no_PH207_model_absent_15x_order <- ggplot(data = PH207_pred_matrix_15x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()

pdf("15x_PH207_reseq_on_no_PH207_model_15x_order.pdf",width=8, height=5)
grid.arrange(PH207_seq_on_no_PH207_model_present_15x_order,PH207_seq_on_no_PH207_model_absent_15x_order,nrow=1)
dev.off()

ref_genome =c("B73v4","Mo17","W22v12")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(PH207_pred_matrix_15x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = PH207_pred_matrix_15x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = PH207_pred_matrix_15x_order %>% filter(absent >=j & Ref_genome ==i)
    total_classifed = as.data.frame(dim(pred_present))[1,] + as.data.frame(dim(pred_absent))[1,]
    present_call = as.data.frame(table(pred_present$TE_PAV))
    false_pos = present_call[1,2]/(present_call[1,2]+present_call[2,2])
    true_pos = present_call[2,2]/(present_call[1,2]+present_call[2,2])
    absent_call= as.data.frame(table(pred_absent$TE_PAV))
    false_neg = absent_call[2,2]/(absent_call[1,2]+absent_call[2,2])
    true_neg = absent_call[1,2]/(absent_call[1,2]+absent_call[2,2])
    pct_captured=total_classifed/total_count 
    QC_stats= as.data.frame(rbind(true_pos,false_pos,true_neg,false_neg,pct_captured))
    QC_stats$Prob_cutoff <- j
    QC_stats$Ref_genome <- i
    colnames(QC_stats[1]) <- "rate"
    QC_stats$metric <- c("true_pos","false_pos","true_neg","false_neg","pct_captured") 
    print(QC_stats)
    out=rbind(out,QC_stats)
  }
}

as.data.frame(out) %>% write.csv("QC_PH207_reseq_on_no_PH207_model_matrix_15x_order.csv")
visualize_QC_PH207_no_PH207_model_15x_order <- read.csv(file ="QC_PH207_reseq_on_no_PH207_model_matrix_15x_order.csv", header= TRUE)
visualize_QC_PH207_no_PH207_model_15x_order$metric <- factor(visualize_QC_PH207_no_PH207_model_15x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))
PH207_reseq_QC_15x_order = visualize_QC_PH207_no_PH207_model_15x_order %>% 
  ggplot(aes(Prob_cutoff,V1,colour = metric))+ facet_wrap(~Ref_genome) + 
  geom_point(position = position_dodge(width = 0.05)) + 
  xlab("Probability Cutoff") +ylab ("Percentage") + 
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#332288","#CC79A7"), 
                     name="Metrics",
                     breaks=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"),
                     labels=c("False Negative", "False Positive","True Negative","True Positive","Total Pecentage of TE Captured")) + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=12),
        legend.text = element_text(size=12),plot.title = element_text(hjust = 0.5)) + ggtitle("PH207 Resequencing Reads") + 
  theme(axis.title.x=element_blank())



######### W22 ############
# using no_W22 training model to predict W22 resequencing reads
predict_prob_W22_15x_order <- predict(W22_training_model_15x_order, newdata=W22_tester_15x, type="prob")
W22_pred_matrix_15x_order <- cbind(W22_tester_15x,predict_prob_W22_15x_order)

W22_seq_on_no_W22_model_present_15x_order <- ggplot(data = W22_pred_matrix_15x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()

W22_seq_on_no_W22_model_absent_15x_order <- ggplot(data = W22_pred_matrix_15x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()
pdf("15x_W22_reseq_on_no_W22_model_15x_order.pdf",width=8, height=5)
grid.arrange(W22_seq_on_no_W22_model_present_15x_order,W22_seq_on_no_W22_model_absent_15x_order,nrow=1)
dev.off()

ref_genome =c("B73v4","Mo17","PH207")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(W22_pred_matrix_15x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = W22_pred_matrix_15x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = W22_pred_matrix_15x_order %>% filter(absent >=j & Ref_genome ==i)
    total_classifed = as.data.frame(dim(pred_present))[1,] + as.data.frame(dim(pred_absent))[1,]
    present_call = as.data.frame(table(pred_present$TE_PAV))
    false_pos = present_call[1,2]/(present_call[1,2]+present_call[2,2])
    true_pos = present_call[2,2]/(present_call[1,2]+present_call[2,2])
    absent_call= as.data.frame(table(pred_absent$TE_PAV))
    false_neg = absent_call[2,2]/(absent_call[1,2]+absent_call[2,2])
    true_neg = absent_call[1,2]/(absent_call[1,2]+absent_call[2,2])
    pct_captured=total_classifed/total_count 
    QC_stats= as.data.frame(rbind(true_pos,false_pos,true_neg,false_neg,pct_captured))
    QC_stats$Prob_cutoff <- j
    QC_stats$Ref_genome <- i
    colnames(QC_stats[1]) <- "rate"
    QC_stats$metric <- c("true_pos","false_pos","true_neg","false_neg","pct_captured") 
    print(QC_stats)
    out=rbind(out,QC_stats)
  }
}

as.data.frame(out) %>% write.csv("QC_W22_reseq_on_no_W22_model_matrix_15x_order.csv")
visualize_QC_W22_no_W22_model_15x_order <- read.csv(file ="QC_W22_reseq_on_no_W22_model_matrix_15x_order.csv", header= TRUE)
visualize_QC_W22_no_W22_model_15x_order$metric <- factor(visualize_QC_W22_no_W22_model_15x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))
W22_reseq_QC_15x_order = visualize_QC_W22_no_W22_model_15x_order %>% 
  ggplot(aes(Prob_cutoff,V1,colour = metric))+ facet_wrap(~Ref_genome) + 
  geom_point(position = position_dodge(width = 0.05)) + 
  xlab("Probability Cutoff") +ylab ("Percentage") + 
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#332288","#CC79A7"), 
                     name="Metrics",
                     breaks=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"),
                     labels=c("False Negative", "False Positive","True Negative","True Positive","Total Pecentage of TE Captured")) + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.title=element_text(size=12),
        legend.text = element_text(size=12),plot.title = element_text(hjust = 0.5)) + ggtitle("W22 Resequencing Reads") + 
  theme(axis.title.x=element_blank())

pdf("order_based_15x_QC.pdf",width=10, height=12)
grid_arrange_shared_legend(B73_reseq_QC_15x_order,Mo17_reseq_QC_15x_order,PH207_reseq_QC_15x_order,W22_reseq_QC_15x_order,nrow=4,ncol = 1) 
dev.off()
setwd("~/Desktop/TE paper/Revision_MS/data_visulization/")
library(tidyverse)
library(lemon)
######### Mo17 ############
# using no_Mo17 training model to predict Mo17 resequencing reads
predict_prob_Mo17_30x_order <- predict(PH207_training_model_30x_all_order, newdata=Mo17_tester_30x, type="prob")
Mo17_pred_matrix_30x_order <- cbind(Mo17_tester_30x,predict_prob_Mo17_30x_order)

# visualize distribution 
Mo17_seq_on_PH207_model_present_30x_order <- ggplot(data = Mo17_pred_matrix_30x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()
Mo17_seq_on_PH207_model_absent_30x_order <- ggplot(data = Mo17_pred_matrix_30x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()
pdf("30x_Mo17_reseq_on_PH207_model_30x_order.pdf",width=8, height=5)
grid.arrange(Mo17_seq_on_PH207_model_present_30x_order,Mo17_seq_on_PH207_model_absent_30x_order,nrow=1)
dev.off()

ref_genome =c("B73v4","PH207","W22v12")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(Mo17_pred_matrix_30x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = Mo17_pred_matrix_30x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = Mo17_pred_matrix_30x_order %>% filter(absent >=j & Ref_genome ==i)
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

as.data.frame(out) %>% write.csv("QC_Mo17_seq_on_PH207_model_matrix_30x_order.csv")
visualize_QC_Mo17_seq_on_PH207_model_30x_order <- read.csv(file ="QC_Mo17_seq_on_PH207_model_matrix_30x_order.csv", header= TRUE)
visualize_QC_Mo17_seq_on_PH207_model_30x_order$metric <- factor(visualize_QC_Mo17_seq_on_PH207_model_30x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))

Mo17_reseq_QC_30x_order = visualize_QC_Mo17_seq_on_PH207_model_30x_order %>% 
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
# using Mo17 training model to predict PH207 resequencing reads
predict_prob_PH207_30x_order <- predict(Mo17_training_model_30x_all_order, newdata=PH207_tester_30x, type="prob")
PH207_pred_matrix_30x_order <- cbind(PH207_tester_30x,predict_prob_PH207_30x_order)

PH207_seq_on_PH207_model_present_30x_order <- ggplot(data = PH207_pred_matrix_30x_order, aes(x = present)) + geom_histogram()+ 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Present") + scale_y_continuous("Number of TE")+theme_bw()

PH207_seq_on_PH207_model_absent_30x_order <- ggplot(data = PH207_pred_matrix_30x_order, aes(x = absent)) + geom_histogram() + 
  facet_wrap(~Ref_genome) + scale_x_continuous(name="Probility of Absent") + scale_y_continuous("Number of TE")+theme_bw()

pdf("30x_PH207_reseq_on_Mo17_model_30x_order.pdf",width=8, height=5)
grid.arrange(PH207_seq_on_PH207_model_present_30x_order,PH207_seq_on_PH207_model_absent_30x_order,nrow=1)
dev.off()

ref_genome =c("B73v4","Mo17","W22v12")
cutoff_group = c("0.6","0.65","0.7","0.75","0.8","0.85","0.9")
out=NULL
for (i in ref_genome){
  for (j in cutoff_group ) {
    total_count = as.data.frame(dim(PH207_pred_matrix_30x_order %>% filter(Ref_genome ==i)))[1,]
    pred_present = PH207_pred_matrix_30x_order %>% filter(present >=j & Ref_genome ==i)
    pred_absent = PH207_pred_matrix_30x_order %>% filter(absent >=j & Ref_genome ==i)
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

as.data.frame(out) %>% write.csv("QC_PH207_seq_on_PH207_modelmatrix_30x_order.csv")
visualize_QC_PH207_seq_on_PH207_model_30x_order <- read.csv(file ="QC_PH207_seq_on_PH207_modelmatrix_30x_order.csv", header= TRUE)
visualize_QC_PH207_seq_on_PH207_model_30x_order$metric <- factor(visualize_QC_PH207_seq_on_PH207_model_30x_order$metric, level=c("false_neg", "false_pos","true_neg","true_pos","pct_captured"))

PH207_reseq_QC_30x_order = visualize_QC_PH207_seq_on_PH207_model_30x_order %>% 
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

pdf("~/Desktop/order_based_30x_QC.pdf",width=10, height=6)
grid_arrange_shared_legend(Mo17_reseq_QC_30x_order,PH207_reseq_QC_30x_order,nrow=2,ncol = 1) 
dev.off()

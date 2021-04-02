library("caret")
library("randomForest")
library("kernlab")
library("tidyverse")
library("caTools")
library("tidyr")
library("tidyverse")

# read in matrix 
reshape_df <- read.csv("30x_matrix_all.csv", header=TRUE)

# from now, only keep TE_id start, and stop, pav information
reshape_df_base = reshape_df %>% select("TE_ID","start","end","TE_order","TE_PAV")
# the whole set will be used for the final training set for all the 30x panel prediction

# downsample training_no_Mo17_ref in order to create a balanced classifier
balanced_reshape_df_base = downSample(reshape_df_base[,2:4], reshape_df_base$TE_PAV, list = FALSE, yname = "TE_PAV")
table(balanced_reshape_df_base$TE_PAV) %>% prop.table()
dim(balanced_reshape_df_base)
# total of 626262 data point are used for the training, subsize it to 500,000 for training

# caret random forest
# 10 folds  3 repeat three times
# take 0.8 of the whole set for training

set.seed(123)
split = sample.split(balanced_reshape_df_base$TE_PAV, SplitRatio = 0.353065 )
training_set = subset(balanced_reshape_df_base, split == TRUE)
write.csv(training_set,file="PH207_Mo17_training_set_mtry123_order.csv")
test_set_within = subset(balanced_reshape_df_base, split == FALSE)


x <- training_set[,1:3]
y <- training_set[,4]

control <- trainControl(method='repeatedcv',
                        number=10,
                        repeats=3)
metric <- "Accuracy"
set.seed(123)
#Number randomely variable selected is mtry
# this will provide information on what is the best mtry for the modeling to get the max accuracy 
tunegrid <- expand.grid(.mtry=c(1:3))
rf_default <- train(TE_PAV~.,
                    data=training_set,
                    method='rf',
                    metric='Accuracy',
                    tuneGrid=tunegrid,
                    trControl=control)
print(rf_default)

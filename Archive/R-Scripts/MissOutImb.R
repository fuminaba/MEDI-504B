# Importing the data an wrangling----

library(tidyverse) # for tidy data analysis
library(readr)     # for fast reading of input files

bc_data0 <-  read.csv(paste0("http://archive.ics.uci.edu/ml/machine-learning-databases/","breast-cancer-wisconsin/breast-cancer-wisconsin.data"), header = FALSE, stringsAsFactors = F)

names (bc_data0) <-  c("sample_code_number", 
                       "clump_thickness", 
                       "uniformity_of_cell_size", 
                       "uniformity_of_cell_shape", 
                       "marginal_adhesion", 
                       "single_epithelial_cell_size", 
                       "bare_nuclei", 
                       "bland_chromatin", 
                       "normal_nucleoli", 
                       "mitosis", 
                       "classes")

str(bc_data0)
bc_data0$bare_nuclei = as.integer(bc_data0$bare_nuclei)
bc_data1 <- bc_data0 %>%
  dplyr::mutate(classes = ifelse(classes == "2", "benign",
                                 ifelse(classes == "4", "malignant", NA)))

# De-duplicate observations----
bc_data2 <- bc_data1 %>% distinct(sample_code_number,.keep_all = TRUE)
row.names(bc_data2) <- bc_data2$sample_code_number
bc_data3 <- bc_data2 %>% select(-sample_code_number)

# Split data into training and testing
set.seed(123)
train.index <- caret::createDataPartition(bc_data3$classes, p = .7, list=FALSE)

train <- bc_data3[ train.index,]
valid  <- bc_data3[-train.index,]

# Check the missing values in the data
summary(train)
sapply(train, function(x) sum(is.na(x)))
dim (train)


MP_plot_train <- VIM::aggr(train, col=c('red','blue'),
                           numbers=TRUE, sortVars=TRUE,
                           labels=names(train), cex.axis=.7,
                           gap=3, 
                           ylab=c("Proportion of Missing Data","Pattern"))

#	Pairwise deletion: ignore missing values
mean(train$bare_nuclei)
mean(train$bare_nuclei, na.rm=TRUE)
mean(train$clump_thickness, na.rm=TRUE)
cor(train$bare_nuclei,train$clump_thickness)
cor(train$bare_nuclei,train$clump_thickness, use = "pairwise.complete.obs")

# Listwise deletion: Remove  the rows where missing values occur. 
listwise <- train[complete.cases(train),]
sapply(listwise, function(x) sum(is.na(x)))
dim(listwise)
mean(listwise$bare_nuclei)
mean(listwise$clump_thickness)
cor(listwise$bare_nuclei,listwise$clump_thickness)

# Impute missing value with mean
# Multiple imputation by chained equations (MICE)

Mean_imp <- mice::mice(train, m=5, method = 'mean', print = FALSE)
mice::densityplot(Mean_imp)
# MULTIPLE IMPUTATION USING PREDICTIVE MEAN MATCHING
pmm_imp <- mice::mice(train,m=10, method = "pmm", print = FALSE)
mice::densityplot(pmm_imp)

imputed<-complete(data=pmm_imp)
mean(imputed$bare_nuclei)
mean(listwise$bare_nuclei)
mean(imputed$clump_thickness)
mean(listwise$clump_thickness)
mean(train$clump_thickness)

# Will not impute validation
sapply(valid, function(x) sum(is.na(x)))
completeValid <- valid[complete.cases(valid),]
table(completeValid$classes)
imputed$classes <- as.factor(imputed$classes)
completeValid$classes <- as.factor(completeValid$classes)

# Class Imbalance
table(imputed$classes)


set.seed(123)
ctrl = trainControl(method = "repeatedcv", 
                                number = 5, 
                                repeats = 3, 
                                summaryFunction = twoClassSummary,
                                classProbs = TRUE,
                                verboseIter = FALSE)

orig_fit <- caret::train(classes ~ .,
                  data = imputed,
                  method = "glm",
                  family = "binomial",
                  preProcess = c("scale", "center"),
                  metric ="ROC",
                  trControl = ctrl)

# Build custom AUC function to extract AUC
# from the caret model object
test_roc <- function(model, data) {
  roc(data$classes,
      predict(model, data, type = "prob")[, "malignant"])
}


pred_original <- predict(orig_fit, completeValid)
cm_orig <- caret::confusionMatrix(pred_original, completeValid$classes, positive="malignant")
result_orig <- c(cm_orig$overall["Accuracy"],cm_orig$byClass['Sensitivity'], cm_orig$byClass['Specificity'], cm_orig$byClass['F1'])
result_orig

test_roc(orig_fit,  completeValid) %>% auc()

# Downsampling
set.seed(123)

ctrl$sampling <- "down"
down_fit <- train(classes ~ .,
                  data = imputed,
                  method = "glm",
                  family = "binomial",
                  preProcess = c("scale", "center"),
                  metric = "ROC",
                  trControl = ctrl)

pred_down <- predict(down_fit, completeValid)
cm_down <- confusionMatrix(pred_down, completeValid$classes, positive="malignant")
result_down <- c(cm_down$overall["Accuracy"],cm_down$byClass['Sensitivity'], cm_down$byClass['Specificity'],  cm_down$byClass['F1'])
result_down


all_results <- data.frame(rbind(result_orig, result_down))
all_results

# Upsampling
set.seed(123)

ctrl$sampling <- "up"
up_fit <- train(classes ~ .,
                  data = imputed,
                  method = "glm",
                  family = "binomial",
                  preProcess = c("scale", "center"),
                  metric = "ROC",
                  trControl = ctrl)

pred_up <- predict(up_fit, completeValid)
cm_up <- confusionMatrix(pred_up, completeValid$classes, positive="malignant")
result_up<- c(cm_up$overall["Accuracy"],cm_up$byClass['Sensitivity'], cm_up$byClass['Specificity'],  cm_up$byClass['F1'])
result_up


all_results <- data.frame(rbind(result_orig, result_down, result_up))
all_results

# SMOTE
set.seed(123)

ctrl$sampling <- "smote"

smote_fit <- train(classes ~ .,
                   data = imputed,
                   method = "glm",
                   family = "binomial",
                   preProcess = c("scale", "center"),
                   metric = "ROC",
                   trControl = ctrl)

pred_smote <- predict(smote_fit, completeValid)
cm_smote <- confusionMatrix(pred_smote, completeValid$classes, positive="malignant")
result_smote <- c(cm_smote$overall["Accuracy"],cm_smote$byClass['Sensitivity'], cm_smote$byClass['Specificity'],  cm_smote$byClass['F1'])
result_smote


all_results <- data.frame(rbind(result_orig, result_down, result_up, result_smote))
all_results

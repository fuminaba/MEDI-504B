# Read clean data
library(tidyverse)
library(caret)
library(pROC)

# Build custom AUC function to extract AUC
# from the caret model object
eval_mod <- function(model, data) {
  pred <- predict(model, data)
  cm <- caret::confusionMatrix(pred, data$classes, positive="malignant")
  auc <- roc(data$classes,
             predict(model, data, type = "prob")[, "malignant"]) %>% auc()
  result <- c(cm$overall["Accuracy"],cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['F1'],AUC=auc)
  return(result)
}


bc_data <- readRDS("bc_clean.RDS")
bc_data$classes <- as.factor(bc_data$classes)
# Read clean data
library(tidyverse)
library(caret)
library(pROC)

# Build custom AUC function to extract AUC
# from the caret model object
eval_mod <- function(model, data) {
  pred <- predict(model, data)
  cm <- caret::confusionMatrix(pred, data$classes, positive="malignant")
  auc <- roc(data$classes,
             predict(model, data, type = "prob")[, "malignant"]) %>% auc()
  result <- c(cm$overall["Accuracy"],cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['F1'],AUC=auc)
  return(result)
}


bc_data <- readRDS("bc_clean.RDS")
bc_data$classes <- as.factor(bc_data$classes)

set.seed(2024)
index <- caret::createDataPartition(bc_data$classes, p = 0.7, list = FALSE)

train_data <- bc_data[index, ]
test_data  <- bc_data[-index, ]

set.seed(859)

cctrl <- trainControl(method = "repeatedcv", 
                      number = 5, 
                      repeats = 3,  
                      savePredictions = TRUE,
                      summaryFunction = twoClassSummary,
                      classProbs = TRUE)
lambda <- 10^seq(-3, 3, length = 100)
alpha <- seq(0, 0.5, length = 6)

set.seed(859)
caret_full <- train(classes ~., 
                    data = train_data, 
                    method = "glm",
                    family="binomial",
                    preProcess = c("center", "scale"),
                    trControl = cctrl)
full <- eval_mod(caret_full, test_data)
full

caret_ridge <- train(classes ~., 
               data = train_data, 
               method = "glmnet",
               family="binomial",
               preProcess = c("center", "scale"),
               trControl = cctrl,
               tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)
caret_ridge$bestTune 

plot(caret_ridge)

names(caret_ridge)
plot(caret_ridge$results$lambda,caret_ridge$results$Sens)

ridge <- eval_mod(caret_ridge, test_data)
ridge
set.seed(859)
caret_lasso <- train(classes ~., 
                     data = train_data, 
                     method = "glmnet",
                     family="binomial",
                     trControl = cctrl,
                     tuneGrid = expand.grid(alpha = 1, lambda = lambda)
)
lasso <- eval_mod(caret_lasso, test_data)
lasso

set.seed(859)
caret_enet <- train(classes ~., 
                     data = train_data, 
                     method = "glmnet",
                     family="binomial",
                     trControl = cctrl,
                     tuneGrid = expand.grid(alpha = alpha, lambda = lambda)
)
enet <- eval_mod(caret_enet, test_data)
enet

rbind(full,ridge, lasso, enet)

# Model coefficients
coefs <- cbind(coef(caret_full$finalModel) %>% as.matrix(),
coef(caret_ridge$finalModel, caret_ridge$bestTune$lambda),
coef(caret_lasso$finalModel, caret_lasso$bestTune$lambda),
coef(caret_enet$finalModel, caret_enet$bestTune$lambda))
colnames(coefs) <- c("full","ridge","lasso","enet")
coefs


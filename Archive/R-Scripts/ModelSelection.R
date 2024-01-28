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

#Try doing backward selection using AIC first, then using BIC. Do you get the same result?
  
M0 = glm(classes ~ 1, data = train_data, family = binomial)  # Null model
M1 = glm(classes ~ ., data = train_data, family= binomial)  # Full model
step.back.aic <- stats::step(M1, direction = "backward", trace = TRUE, k = 2)
step.back.bic <-  stats::step(M1, direction = "backward", trace = TRUE, k = log(dim(train_data)[1]))
summary(step.back.bic)

#Now Forward
step.fwd.aic <-  stats::step(M0, scope = list(lower = M0, upper = M1), direction = "forward", trace = FALSE, k = 2)
summary(step.fwd.aic)
step.fwd.bic <-  stats::step(M0,scope = list(lower = M0, upper = M1), direction = "forward", trace = FALSE, k = log(dim(train_data)[1]))
summary(step.fwd.bic)
## Exhaustive searches
library(bestglm)
x_tr <- train_data 
x_tr$classes <- ifelse(train_data$classes=="malignant",1,0)

res.bestglm <-bestglm::bestglm(Xy = x_tr,
                               family = binomial,
                               IC = "AIC",     # Information criteria for
                               method = "exhaustive")

## Show top 5 models
res.bestglm$BestModels

summary(res.bestglm$BestModel)

# in caret

trainX <- train_data[,-10]
trainY <- train_data$classes

cctrl <- trainControl(method = "repeatedcv", 
                       number = 5, 
                       repeats = 3,  
                       savePredictions = TRUE,
                       summaryFunction = twoClassSummary,
                      classProbs = TRUE
)
set.seed(849)
caret_full <- train(classes~.,
                   data = train_data,
                   method = "glm", 
                   family = "binomial",
                   trControl = cctrl,
                   preProc = c("center", "scale"),
                   metric = "ROC",
                   trace = 0)
full <- eval_mod(caret_full, test_data)

caret_fwd <- train(classes~.,
                             data = train_data,
                             method = "glmStepAIC", 
                             direction = "forward",
                             trControl = cctrl,
                             preProc = c("center", "scale"),
                             metric = "ROC",
                             trace = 0)
fwd <- eval_mod(caret_fwd, test_data)

caret_back <- train(classes~.,
                   data = train_data,
                   method = "glmStepAIC", 
                   direction = "backward",
                   trControl = cctrl,
                   preProc = c("center", "scale"),
                   metric = "ROC",
                   trace = 0)
back <- eval_mod(caret_back, test_data)
?MASS::stepAIC
?stats::step

caret_both <- train(classes~.,
                    data = train_data,
                    method = "glmStepAIC", 
                    direction = "both",
                    trControl = cctrl,
                    preProc = c("center", "scale"),
                    metric = "ROC",
                    trace = 0)
both <- eval_mod(caret_both, test_data)
rbind(full, fwd, back, both)

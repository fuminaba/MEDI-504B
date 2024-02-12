# Read clean data
library(tidyverse)
library(caret)
library(pROC)

rootfolder <- 'F:/Courses/MEDI_504B/springt_2/'

# Load data
diabetic_data0 <- 
  read.csv("F:/Courses/MEDI_504B/springt_2/Datasets/diabetic_retinopathyDataSet_train.csv")

# >>> Assign feature names <<< #
names(diabetic_data0) <-  
  c("quality", "pre_screening", 
    "ma1", "ma2", "ma3", "ma4", "ma5", "ma6", 
    "exudate1", "exudate2", "exudate3", "exudate4",
    "exudate5", "exudate6", "exudate7","exudate8",
    "macula_opticdisc_distance", "opticdisc_diameter",
    "am_fm_classification", "classes")

# >>> Recode numeric outcome variable to string, then factor <<< #
diabetic_data1 <- diabetic_data0 %>%
  mutate(classes = ifelse(classes == 0, "No",
                          ifelse(classes == 1, "Sign", NA))) %>%
  mutate(classes = as_factor(classes))

# Remove bad image quality data
diabetic_data1 <- diabetic_data1[complete.cases(diabetic_data1),] %>% 
  filter(quality == 1) %>% 
  select(-quality)

# Remove "quality" as a variable
diabetic_data1 <- diabetic_data1[,-1]

# Build custom AUC function to extract AUC
# from the caret model object
eval_mod <- function(model, data) {
  pred <- predict(model, data)
  cm <- caret::confusionMatrix(pred, data$classes, positive="Sign")
  auc <- roc(data$classes,
             predict(model, data, type = "prob")[, "Sign"]) %>% auc()
  result <- c(cm$overall["Accuracy"],cm$byClass['Sensitivity'], cm$byClass['Specificity'], cm$byClass['F1'],AUC=auc)
  return(result)
}

# Split training and testing data
set.seed(2024)
index <- caret::createDataPartition(diabetic_data1$classes, p = 0.7, list = FALSE)
train_data <- diabetic_data1[index, ]
test_data  <- diabetic_data1[-index, ]
train_data$classes %>% table(.)
set.seed(2024)

# Decision Tree
library(rpart)
tree_mod <- rpart(
  formula = classes ~. ,
  data    = train_data,
  method  = "class"
)
tree_mod
library(rpart.plot)
rpart.plot(tree_mod)


# Pruning
# We prune the tree using the cost complexity criterion. Can a less deep tree give comparable results. If so, go with the shallower tree to reduce the likelihood of overfitting. We used the a parameter called complexity parameter (CP)

printcp(tree_mod)

# it computes the cross validation error for each value of cp CP=0.01 gives the lowest error
opt <- which.min(tree_mod$cptable[,"xerror"])
#get its value
cp <- tree_mod$cptable[opt, "CP"]

# we can prune the tree based on this CP
pruned_model <- prune(tree_mod,cp)
#plot tree
rpart.plot(pruned_model)
#Note that rpart will use a default CP value of 0.01 if you don’t specify one in prune.

pred <- predict(object = tree_mod,   # model object 
                newdata = test_data,
                type="class")  # test dataset
caret::confusionMatrix(pred,test_data$classes)

set.seed(2024)
caret_tree <- train(
  classes ~ .,
  data = train_data,
  method = "rpart",
  metric ="Sens",
  trControl = trainControl(method = "cv", number = 5,
                           classProbs = T, summaryFunction = twoClassSummary),
  tuneLength = 20
)
ggplot(caret_tree)

caret_tree$bestTune
rpart.plot(caret_tree$finalModel)

tree <- eval_mod(caret_tree,test_data)

# Treebag
tree_bag1 <- ipred::bagging(
formula = classes ~ .,
data = train_data,
nbagg = 500,  
coob = TRUE
)

tree_bag1
set.seed(2024)
caret_bag <- train(
  classes ~ .,
  data = train_data,
  method = "treebag",
  trControl = trainControl(method = "cv", number = 5,
                           classProbs = T, summaryFunction = twoClassSummary),
  metric ="Sens",
  nbagg = 20,  
  control = rpart.control(minsplit = 2, cp = 0)
)
caret_bag
bag <- eval_mod(caret_bag,test_data)

rbind(tree, bag)

# Random Forest
library(randomForest)
# Train a Random Forest
set.seed(2024)  # for reproducibility
rf_model <- randomForest(formula = classes ~ ., 
                         data = train_data)

# Print the model output                             
print(rf_model)

caret_rf <- train(
  classes ~ .,
  data = train_data,                         
  method = "ranger",
  metric = "Sens",
  trControl = trainControl(method = "cv", number = 5,
                           classProbs = TRUE, summaryFunction = twoClassSummary),
  importance = "impurity"
)

forest <- eval_mod(caret_rf,test_data)
rbind(tree, bag, forest)

# Using grid search
control <- trainControl(method="cv", number=5, search="grid")
set.seed(2024)
tunegrid <- expand.grid(.mtry=c(5:15))

train_ctrl <- trainControl(method="cv", 
                           number=5, 
                           search = "grid" ,
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary,
                           selectionFunction = "best",
                           returnResamp = "all",
                           savePredictions = TRUE,
                           verboseIter = TRUE
)

model_cv_grid <- train(classes ~ .,
                       data = train_data,
                       method = "rf", 
                       metric = "Sens", 
                       trControl = train_ctrl, 
                       tuneGrid = tunegrid
) 

model_cv_grid
plot(model_cv_grid)

optimal_mtry <- model_cv_grid$bestTune$mtry
optimal_rf <- randomForest(classes ~ .,
                           data = train_data,
                           mtry = optimal_mtry,
                           importance = TRUE, 
                           keep.forest = TRUE
)
rf_grid <- eval_mod(optimal_rf, test_data)
rbind(tree, bag, forest, rf_grid)

varImpPlot(optimal_rf)


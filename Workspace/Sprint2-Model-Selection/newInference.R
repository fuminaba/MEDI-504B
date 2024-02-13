# =============================== #
# >>> Define Useful Functions <<< #
# =============================== #
library(readr)
library(here)
library(bestglm)
library(caret)
library(pROC)
library(dplyr)
library(tidyverse)
# >>> Data loading function <<< #
data_load <- function(data.path) {
    
    # This function loads and processes the data for inference of 
    # tree-based models. 
    #   It takes in a string (data.path) which is passed to read.csv. 
    #   It assigns feature names to columns, converts output variable to a 
    #   factor. Then converts binary features to factors, and removes the 
    #   'quality' feature.

    data.table <- data.path %>%
        read.csv(header = TRUE, stringsAsFactors = F)
    names(data.table) <-  
        c("quality", "pre_screening", 
          "ma1", "ma2", "ma3", "ma4", "ma5", "ma6", 
          "exudate1", "exudate2", "exudate3", "exudate4",
          "exudate5", "exudate6", "exudate7","exudate8",
          "macula_opticdisc_distance", "opticdisc_diameter",
          "am_fm_classification", "classes")
    
    data.table <- data.table %>%
        mutate(classes = ifelse(classes == 0, "No",
                                ifelse(classes == 1, "Sign", NA))) %>%
        mutate(classes = as_factor(classes),
               pre_screening = as_factor(pre_screening)) %>%
        filter(quality == 1) %>%
        select(-quality)
    return (data.table)
}

# >>> Model Evaluation <<< #
eval_mod <- function(model, data) {
    # Function by Dr.Aline Talhouk to compute a number of evaluation metrics.
    #   Computes the following: 
    #       * Accuracy
    #       * Sensitivity (metric of interest)
    #       * Specificity 
    #       * F1 score
    pred <- predict(model, data)
    cm <- caret::confusionMatrix(pred, data$classes, positive="Sign")
    print(cm)
    auc <- pROC::roc(data$classes,
                     predict(model, data, type = "prob")[, "Sign"]) %>% 
        auc()
    result <- c(cm$overall["Accuracy"],
                cm$byClass['Sensitivity'], 
                cm$byClass['Specificity'], 
                cm$byClass['F1'],AUC=auc)
    return(result)
}

# ================================= #
# >>> Load necessary components <<< #
# ================================= #

# >>> Load data <<< #
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # >>> This may need to be changed to correct path
new.data.path <- '../../Datasets/diabetic_retinopathyDataSet_train.csv' # >>> Change this to path to test data 
test.data <- data_load(new.data.path)

# >>> Load model(s) <<< #
xgb_tree <- readRDS("./xgboost_tree.RDS")

# >>> Evaluate the model <<< #
xgb_eval <- eval_mod(xgb_tree, test.data)
xgb_eval


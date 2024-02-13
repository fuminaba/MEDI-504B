data_load <- function(data.path) {
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

data_preproc <- function(datatable) {
    datatable <- datatable %>%
        dplyr::filter(quality == 1) %>%
        dplyr::mutate(pre_screening = as_factor(pre_screening),
                      am_fm_classification = as_factor(am_fm_classification)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(ma_mean = (mean(ma1, ma2, ma3, ma4, ma5, ma6)),
                      exudate1 = log(exudate1 + 1),
                      exudate_mean = mean(exudate4, exudate5, exudate6, exudate7, exudate8) + 1,
        ) %>% 
        dplyr::select(-c(quality, exudate2, exudate3, exudate4,
                         exudate5, exudate6, exudate7, exudate8,
                         ma1, ma2, ma3, ma4, ma5, ma6,
                         macula_opticdisc_distance, opticdisc_diameter,
                         am_fm_classification))
    
    return(datatable)
}

eval_mod <- function(model, data) {
    # Function by Dr.Aline Talhouk
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

linearity_plot <- function(model, data) {
    preds <- predict(model, data, type = 'prob')
    logodds_outcome <- log(preds / (1-preds))
    return (logodds_outcome)
}
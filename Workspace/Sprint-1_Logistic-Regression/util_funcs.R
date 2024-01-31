data_preproc <- function(datatable) {
    datatable <- datatable %>%
        dplyr::filter(quality == 1) %>%
        dplyr::mutate(pre_screening = as_factor(pre_screening),
                      am_fm_classification = as_factor(am_fm_classification)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(ma_sum = mean(ma1, ma2, ma3, ma4, ma5, ma6),
                      exudate_mean = mean(exudate5, exudate6, exudate7, exudate8)) %>% 
        dplyr::select(-c(quality, exudate1, exudate2, exudate3, 
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
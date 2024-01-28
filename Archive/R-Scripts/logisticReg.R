# Read clean data
library(tidyverse)

bc_data <- read.csv(paste0("http://archive.ics.uci.edu/ml/machine-learning-databases/","breast-cancer-wisconsin/breast-cancer-wisconsin.data"), header = FALSE, stringsAsFactors = F)
names (bc_data) <-  c("sample_code_number", 
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

str(bc_data)
bc_data <- bc_data %>% distinct(sample_code_number,.keep_all = TRUE)



bc_data$bare_nuclei = as.integer(bc_data$bare_nuclei)
# NA introduced by coercion is just a warning message that lets us know that some values got coerced into NA

bc_data <- bc_data %>%
  dplyr::mutate(classes = ifelse(classes == "2", "benign",
                                 ifelse(classes == "4", "malignant", NA)))

str(bc_data)

bc_data <- na.omit(bc_data)
# notice that you get an error requiring the Y values to be between zero and one


# we also don't really care about sample code number, as it is not a biological variable.
row.names(bc_data) <- bc_data$sample_code_number
bc_data$classes <- as.factor(bc_data$classes)

bc_data <- bc_data %>%select(-sample_code_number)
# we reset classes from character type to factor
model_logit <-  glm(classes~., data= bc_data, family = binomial)
summary(model_logit)
confint(model_logit)
# notice that these are on the log odds scale.
# Publication ready results
fit.apparent <- Publish::publish(model_logit, digits=1)$regressionTable

# Assessing the model's accuracy
# Apparent accuracy----
pred_class <- predict(model_logit, bc_data, type = "response")

fitted_mod <- fitted(model_logit)
plot(pred_class, fitted_mod)
abline(0,1, col="red")

# fitted values are the same as predicted values using the training data

ggplot(data.frame(y_linear = pred_class), aes(pred_class)) +
  geom_density() +
  ggtitle("Probabilities predicted by linear model of binary default")+
  xlab("Predicted class")

pred_class <- ifelse(pred_class>0.5,1,0) %>% factor(labels = c("benign","malignant"))
table(pred_class)

## ROC Plots
library(recipes)
library(caret)
library(pROC)
library(plotROC)

caret:: confusionMatrix(pred_class, bc_data$classes, positive ="malignant")
?caret::confusionMatrix
roc_mfit <- roc(bc_data$classes, as.numeric(pred_class))
?roc
# the first argument is the observed response, the second argument is the predicted response. Be careful this is not a norm and you have to check for every function.

plot(roc_mfit, col = "red", print.auc = T)

# Resampling
set.seed(2024)
index <- caret::createDataPartition(bc_data$classes, p = 0.7, list = FALSE)
#By default, createDataPartition does a stratified random split of the data.

train_data <- bc_data[index, ]
test_data  <- bc_data[-index, ]

bind_rows(data.frame(group = "train", train_data),
          data.frame(group = "test", test_data)) %>%
  gather(x, y, clump_thickness:mitosis) %>%
  ggplot(aes(x = y, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)
# can see that the training set and the test said are pretty similar in terms of covariates

set.seed(2024)

logit_cv <- caret::train(classes ~ .,
                         data = train_data,
                         method = "glm",
                         family = "binomial",
                         preProcess = c("scale", "center"),
                         ## Center and scale the predictors for the training
                         ## set and all future samples.
                         trControl = trainControl(method = 'none')
)

summary(logit_cv)

trControl_params = trainControl(method = "repeatedcv", 
                         number = 5, 
                         repeats = 3, 
                         savePredictions = TRUE, 
                         summaryFunction = twoClassSummary,
                         verboseIter = FALSE)

#method specifies how resampling will be done. Examples include cv, boot, LOOCV, repeatedcv, and oob.
#number specifies the number of times resampling should be done for methods that require resample, such as, cv and boot.
#repeats specifies the number of times to repeat resampling for methods such as repeatedcv
?trainControl

logit_cv <- caret::train(classes ~ .,
                         data = train_data,
                         method = "glm",
                         family = "binomial",
                         preProcess = c("scale", "center"),
                         ## Center and scale the predictors for the training
                         ## set and all future samples.
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 5, 
                                                  repeats = 3, 
                                                  savePredictions = TRUE)
)


logit_cv

# we used 441 samples that had binary class benign malignant used tenfold cross validation repeated 10 times. It also provides the cross-validated results accuracy, and Kappa.

names(logit_cv)

#There are many parameters that are stored in the training object. 
logit_cv$results
logit_cv$finalModel

#The finalModel is a model object, in this case, the object returned from glm(). This final model, is fit to all of the supplied training data. This model object is often used when we call certain relevant functions on the object returned by train(), such as summary()

summary(logit_cv$finalModel)
summary(glm(classes~., data= train_data, family = "binomial"))

#To predict new samples, predict can be used. For classification models, the default behavior is to calculate the predicted class. The option type = "prob" can be used to compute class probabilities from the model. For example:

predictions <- predict(logit_cv, newdata = test_data)
str(predictions)
# these predictions are made on the held out data set
caret::confusionMatrix(predictions, test_data$classes)
# first predictions then the reference
?caret::confusionMatrix
caret::confusionMatrix(predictions, test_data$classes, positive = "malignant")


predictions <- predict(logit_cv, newdata = test_data, type="prob")
str(predictions)
rowSums(predictions)
head(predictions)

# what do you notice between the accuracy measures on the test set versus the training set?

results <- data.frame(actual = test_data$classes,
                      predict(logit_cv, test_data, type = "prob"))

results$prediction <- ifelse(results$benign > 0.95, "benign",
                             ifelse(results$malignant > 0.95, "malignant", NA))

results$correct <- ifelse(results$actual == results$prediction, TRUE, FALSE)

ggplot(results, aes(x = prediction, fill = correct)) +
  geom_bar(position = "dodge")


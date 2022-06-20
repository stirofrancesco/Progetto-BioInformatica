library(data.table)
library(devtools)
library(caTools)
library(e1071)
library(caret)

#Caricamento della matrice di conteggi
count <- read.delim("count.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")

#Rimozione delle colonne costanti
count[c("X101011101", "X111101101","X111001111")] <- list(NULL)

#Preparazione della matrice da passare a SVM.
{
  count$SUBTYPE = factor(count$SUBTYPE, levels = c("LumA","Her2","LumB","Normal","Basal"))
  count <- count %>%
    select(SUBTYPE, everything())

  count[2:502]  <- lapply(count[2:502], as.numeric)
}

#Costruzione di training set e test set
set.seed(123)
split = sample.split(count, SplitRatio = 0.8)

training_set = subset(count, split == TRUE)
test_set = subset(count, split == FALSE)


#SVM con metodo Radial Sigma
set.seed(123)
fitc = trainControl(method = "cv", number = 10, search = "random", savePredictions = T)
svm_rad = train(SUBTYPE~ . ,data = training_set, method = "svmRadialSigma",
                trControl=fitc,
                tuneLength = 400,
                preProcess=c("center","scale"))

#Parametro migliore
svm_rad$bestTune

#Accuracy con parametro migliore
svm_rad$results[svm_rad$results$Accuracy == max(svm_rad$results$Accuracy), ]

#Predizione usando il test set
pred <- predict(svm_rad,test_set)

#Matrice di confusione
confusionMatrix(table(test_set[,"SUBTYPE"],pred))


#GRID search basata sulla ricerca random
#set.seed(123)
#fitc = trainControl(method = "cv", number = 5, search = "random", savePredictions = T)
#svm_rad_grid = train(SUBTYPE~ . ,data = count, method = "svmRadialSigma",
#                trControl=fitc,
#                tuneGrid = expand.grid(
#                    .sigma = seq(0.01558077,0.2,length=40),
#                    .C = seq(5,20, length=40)
#                ),preprocess=c("center","scale"))

#best parameter
#svm_rad_grid$bestTune

#accuracy with best parameter
#svm_rad_grid$results[svm_rad_grid$results$Accuracy == max(svm_rad_grid$results$Accuracy), ]

#Confusion matrix with optimal parameters
#sub_svm=subset(svm_rad_grid$pred,(
#  svm_rad_grid$pred$C==svm_rad_grid$bestTune$C &
#    svm_rad_grid$pred$sigma==svm_rad_grid$bestTune$sigma))

#confusionMatrix(table(sub_svm$pred,sub_svm$obs),positive = "1")



#---------------


classifier = svm(formula = SUBTYPE ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'poly')

y_pred = predict(classifier, newdata = test_set)

#Confusion Matrix

expected_value <- factor(test_set$SUBTYPE)

example <- confusionMatrix(data=y_pred, reference = expected_value)
print(example)

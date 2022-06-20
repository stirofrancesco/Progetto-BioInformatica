library(data.table)
library(devtools)
library(caTools)
library(e1071)
library(caret)

count <- read.delim("count.txt",
                    header = TRUE,
                    sep = "\t",
                    colClasses = "character")

count[c("X101011101", "X111101101","X111001111")] <- list(NULL)

count$SUBTYPE = factor(count$SUBTYPE, levels = c("LumA","Her2","LumB","Normal","Basal"))
count <- count %>%
  select(SUBTYPE, everything())

count[2:502]  <- lapply(count[2:502], as.numeric)

set.seed(123)
split = sample.split(count, SplitRatio = 0.80)

training_set = subset(count, split == TRUE)
test_set = subset(count, split == FALSE)

set.seed(123)
fitc = trainControl(method = "cv", number = 10, search = "random", savePredictions = T)
svm_rad = train(SUBTYPE~ . ,data = training_set, method = "svmLinear",
                trControl=fitc,
                preProcess=c("center","scale"))

pred <- predict(svm_rad,test_set)

confusionMatrix(table(test_set[,"SUBTYPE"],pred))



#---------------


classifier = svm(formula = SUBTYPE ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'linear')

y_pred = predict(classifier, newdata = test_set)

#Confusion Matrix

expected_value <- factor(test_set$SUBTYPE)

example <- confusionMatrix(data=y_pred, reference = expected_value)
print(example)
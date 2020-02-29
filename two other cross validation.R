# load the dataset, change all the continuous predictors to numeric, taking out
# missing value, get BreastCancer3 for downstream analysis. last class column should be 
# be a factor
library(mlbench)
library(MASS)
library(caret)
library(e1071)

set.seed(5)
data(BreastCancer)
BreastCancer1 <- BreastCancer[,-1]
BreastCancer_predictor <- apply(BreastCancer1[,-10],2,as.numeric)
BreastCancer_predictor_df <- as.data.frame(BreastCancer_predictor)
BreastCancer2 <- cbind(BreastCancer_predictor_df,BreastCancer1$Class)
colnames(BreastCancer2)[10] <- "class"
row_has_na <- apply(BreastCancer1,1,function(x){any(is.na(x))})
sum(row_has_na)
BreastCancer3 <- BreastCancer2[!row_has_na,]

BreastCancer3$class <- as.factor(BreastCancer3$class)

# using caret package to achieve kfold and leave one out cross validation
kfold <- function(nfold,dataset=BreastCancer3){
  CrossVal <- trainControl(method = "repeatedcv",number=10,savePredictions = TRUE)
  Fit1 <- train(class~.,data=dataset,method="glm",family="binomial",trControl= CrossVal)
  Fit1$results
}

kfold_list <- c(10,15,20,25)
kfold_result <- lapply(kfold_list,kfold)
names(kfold_result) <- c("kfold is 10","kfold is 15","kfold is 20","kfold is 25")


leave_one_out <- function(dataset=BreastCancer3){
  TrainControl <- trainControl(method='LOOCV')
  Fit2 <- train(class~.,data=dataset,trControl=TrainControl,method="glm",family=binomial)
  Fit2$results
}

leave_one_out()


# some useful R code from recent coding:

# scan funciton
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

# perform t-test to find differential expressed gene between two groups
training_t_stat<-data.frame(sapply(seq(nrow(trainingA)), function(x) { abs(as.numeric(t.test(trainingA[x,], trainingB[x,])$statistic)) }))


# from the scratch way to achieve kfold cross validation
Positive_groups <- split(sample(colnames(Positive)), 1+(seq_along(colnames(Positive)) %% nfold))
Negative_groups <- split(sample(colnames(Negative)), 1+(seq_along(colnames(Negative)) %% nfold))

# using glmnet package to perform lasso regression
training <- transpose4logistic(trainingA,trainingB,selected_genes)
testing <- transpose4logistic(testA,testB,selected_genes)
    
    x <- as.matrix(subset(training,select=-c(ER_status)))
    x_pre <- as.matrix(subset(testing,select=-c(ER_status)))
    glmmod <- glmnet(x, y=as.factor(training$ER_status), alpha=1, family="binomial")
    #cv.glmmod <- cv.glmnet(x, y=as.factor(training$ER_status), alpha=1)
    prediction <- round(predict(glmmod, newx = x_pre, s = 0.1),3)
    hard_pred <- ifelse(prediction >= 0.5, "Positive","negative")
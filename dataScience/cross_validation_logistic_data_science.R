# This Rscript contains logistic regression, lasso regression, cross-validation.

# laod data, clean the data for logistic regression
library(mlbench)
library(MASS)
data(BreastCancer)
BreastCancer1 <- BreastCancer[,-1]
BreastCancer_predictor <- apply(BreastCancer1[,-10],2,as.numeric)
BreastCancer_predictor_df <- as.data.frame(BreastCancer_predictor)
BreastCancer2 <- cbind(BreastCancer_predictor_df,BreastCancer1$Class)
colnames(BreastCancer2)[10] <- "class"
row_has_na <- apply(BreastCancer1,1,function(x){any(is.na(x))})
sum(row_has_na)
BreastCancer3 <- BreastCancer2[!row_has_na,]

# ramdom sampling to get training and testing data, not n-fold
"
three way of doing validation:
1. fixed proportion, like 80%/20%
2. k-fold cross validation
3. leave-one-out cross validation
"
set.seed(1)
mis_rate_vector <- c()
for (i in 1:100){
  index <- sample(1:683,546,replace = F)
  training_set <- BreastCancer3[index,]
  testing_set <- BreastCancer3[-index,]
  model <- glm(class~., data=training_set,family = binomial)
  prediction <- round(predict(model,newdata=testing_set,type="response"),3)
  hard_pred <- ifelse(prediction >= 0.5, "malignant","benign")
  misclassification <- table(testing_set$class,hard_pred)
  mis_rate <- 1-sum(diag(misclassification))/137
  mis_rate_vector <- c(mis_rate_vector,mis_rate)
}

# calculate overall performance

norm_confit <- function(vector_input){
  me = mean(vector_input)
  s = sd(vector_input)
  n = length(vector_input)
  se = s/sqrt(n)
  lower_bound = me - 2*se   # se or s all make sense i think, by definition, should use s to get 95% area under the curve, but se would make precise range when having larger sample size, see gmail bookmark
  upper_bound = me + 2*se
  c = c("mean" = me, "sd" = s, "lower_bound" = lower_bound,"upper_bound" = upper_bound)
  return(c)
}

norm_confit(mis_rate_vector)



# a simple logistic regression example, contains some other useful code including model selection function
library(MASS)
data(birthwt)
dim(birthwt)
head(birthwt,6)
?birthwt


birthwt1 <- birthwt[,-10]
class <- apply(birthwt1,2,class)

birthwt1$race <- factor(birthwt1$race,labels=c("white","black","other"))
# this factor function is specially designed for regression, the labels just tell the computer to set "white" as base.
table(birthwt1$race)

model <- glm(low~.,family = binomial,data = birthwt1)
summary(model)

pchisq(201.38,179,lower.tail=F)   # null hypothesis: model fits well, if p>0.05, accept null hypothesis
selection <- stepAIC(model,direction="backward")  # AIC, the less, the better, same for BIC
confint(model)
oddsratio <- exp(model$coefficients)
pred <- predict(model,newdata=birthwt1,type = "response")
hard_pred <- ifelse(pred>=0.5,1,0)
ConfusionMatrix <- table(birthwt1$low,hard_pred)
rownames(ConfusionMatrix) <-c("observed normal","observed low weight")
colnames(ConfusionMatrix) <- c("Predicted normal","predicted low weight")
misclassification <- (36+13)/(117+23+36+13)


# using caret package to perform k-fold cross-validation or leave-one-out cross validation
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


# some useful R code from Project1:

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
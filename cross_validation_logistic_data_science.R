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
  error = qnorm(0.975) * s/sqrt(n)
  lower_bound = me - error
  upper_bound = me + error
  c = c("mean" = me, "sd" = s, "lower_bound" = lower_bound,"upper_bound" = upper_bound)
  return(c)
}

norm_confit(mis_rate_vector)


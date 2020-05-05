"
This Rscript contains how to do multivariate linear regression, how to do classification tree and regression tree
using rpart packages.
"



library(faraway)
library(MASS) # for stepAIC function
library(rpart)
library(rpart.plot)
data(fat)
dim(fat)
head(fat,10)
sapply(fat,class)  # check the type of each predictor, all numeric
summary(fat)
row_na <- apply(fat,1,function(x){any(is.na(x))})
sum(row_na)
fat_2 <- fat[,-c(2,3)]

# multiple regression model, multivariate linear regression
fit <- lm(brozek ~.,data=fat_2 )
summary(fit)
coefficients(fit) # model coefficients
confint(fit, level=0.95) # CIs for model parameters 
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table 
vcov(fit) # covariance matrix for model parameters 
influence(fit) # regression diagnostics

# please refer to this webpage: https://www.statmethods.net/stats/regression.html

# check goodness-of-fit (same as regression tree, using R^2)
## RSS <- cal_RSS(fat_2$brozek,fitted(fit))
## TSS <- cal_TSS(fat_2$brozek)
## R2 <- R_square(RSS,TSS)
## R2

# feature selection backward elimination
selection <- stepAIC(fit,direction="backward")
selection$anova  # display the result
# eliminate one predictor once, see if AIC, negatively correlated with likelihood, decrease. Until AIC won't decrease,
# that will be the final model

# Fit the best model
fit_best <- lm(brozek ~ weight + adipos + free + chest + abdom + thigh + ankle + 
                 biceps + forearm + wrist, data=fat_2)

# build regression tree
MB <- rpart(brozek ~ . , data = fat_2,na.action=na.omit)  # regression tree aims to minimize variance of each leafnode, 
# this is the criteria to split up one node, intead of using GINI or entropy in classification tree.
# rpart will automatically recognize if it is a classification tree or regression tree by detecting brozek's value(discrete or continuous)
# the advantage of it over linear regression is that it offer stratification and more interpretability
rpart.plot(MB, type = 4, extra = 1, col = "red")  # extra control the aesthetic features, mean was reported in each node
MB_pre <- predict(MB,newdata=fat_2)

# performance of the tree (R square  = 1 - RSS/TSS)

# define functions to calculate RSS and TSS
cal_RSS <- function(observed,predicted){
  sum = 0
  for (i in length(observed)){
    diff = (observed[i] - predicted[i])^2
    sum = sum + diff
  }
  return(sum)  
}

cal_TSS <- function(observed){    # TSS = variance * N, TSS = sd^2 * N
  me = mean(observed)
  TSS_array = sapply(observed,function(x){(x-me)^2})
  TSS = sum(TSS_array)
  return(TSS)
}


# define a function to calculate R square
R_square <- function(RSS,TSS){
  result = 1 - RSS/TSS
  return(result)
}

# start to calculate
RSS <- cal_RSS(fat_2$brozek,MB_pre)
TSS <- cal_TSS(fat_2$brozek)
R2 <- R_square(RSS,TSS)

# plot the result
plot(fat_2$brozek, MB_pre, pch = 16, col = "red", xlab = "Observed Fat",
     ylab = "Predicted Fat", main = "Regression Tree Output on the fat data")








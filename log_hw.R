library(MASS)
data(birthwt)
#setwd("C:\\Users\\ligk2E\\Desktop\\coding\\logistic_regression")
dim(birthwt)
head(birthwt,6)

birthwt1 <- birthwt[,-10]
class <- apply(birthwt1,2,class)

birthwt1$race <- as.factor(birthwt1$race)
table(birthwt1$race)

model <- glm(low~.,family = binomial,data = birthwt1)
summary(model)



Data <- data.frame(low=c(0,1),age=c(56,89),lwt=c(89,144),race=as.factor(c(3,2)),smoke=c(0,1),
                   ptl=c(1,0),ht=c(0,0),ui=c(0,0),ftv=c(0,1))

Prediction1 <- predict(model,newdata = Data,type="response")








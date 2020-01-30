library(MASS)
data(birthwt)
dim(birthwt)
head(birthwt,6)
?birthwt

birthwt1 <- birthwt[,-10]
class <- apply(birthwt1,2,class)

birthwt1$race <- factor(birthwt1$race,labels=c("white","black","other"))
table(birthwt1$race)

model <- glm(low~.,family = binomial,data = birthwt1)
summary(model)

pchisq(201.38,179,lower.tail=F)

confint(model)
oddsratio <- exp(model$coefficients)
pred <- predict(model,newdata=birthwt1,type = "response")
hard_pred <- ifelse(pred>=0.5,1,0)
ConfusionMatrix <- table(birthwt1$low,hard_pred)
rownames(ConfusionMatrix) <-c("observed normal","observed low weight")
colnames(ConfusionMatrix) <- c("Predicted normal","predicted low weight")
misclassification <- (36+13)/(117+23+36+13)







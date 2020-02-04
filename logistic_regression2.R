library(mlbench)
data(BreastCancer)
dim(BreastCancer)
head(BreastCancer,6)

BreastCancer1 <- BreastCancer[,-1]
class(BreastCancer1$Cl.thickness)
apply(BreastCancer1,2,class)
typeof(BreastCancer1$Cl.thickness)
mode(BreastCancer1$Cl.thickness)
str(BreastCancer1$Cl.thickness)
sapply(BreastCancer1, class)
lapply(BreastCancer1,class)

sapply(BreastCancer1, class)
lapply(BreastCancer1,class)

BreastCancer_predictor <- apply(BreastCancer1[,-10],2,as.numeric)
BreastCancer_predictor_df <- as.data.frame(BreastCancer_predictor)
BreastCancer2 <- cbind(BreastCancer_predictor_df,BreastCancer1$Class)
colnames(BreastCancer2)[10] <- "class"
row_has_na <- apply(BreastCancer1,1,function(x){any(is.na(x))})
sum(row_has_na)
BreastCancer3 <- BreastCancer2[!row_has_na,]

summary(BreastCancer3)
model <- glm(class~.,family=binomial,data=BreastCancer3)
summary(model)
selection <- stepAIC(model,direction="backward")

confint(model)
oddsratio<-exp(model$coefficients)
pred<-predict(model,newdata=BreastCancer3,type = "response")
hard_pred <- ifelse(pred>0.5,"malignant","benign")
ConfusionMatrix <- table(BreastCancer3$class,hard_pred)
rownames(ConfusionMatrix) <- c("observed benign","observed malignant")
colnames(ConfusionMatrix) <- c("predicted benign","predicted maligant")
misclassification <- 

rm(list=ls())

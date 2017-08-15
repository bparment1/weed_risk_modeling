#MAJ=major invader=1 #MIN=minor invader=0
#68 instances of 1 #68 instances of 0
setwd("/Users/chinchuharris/modeling_weed_risk/data")
data=read.csv(file="secondmodel.csv") #data for prediction
library(randomForest) #for random forests
library(lattice)
library(ggplot2)
library(caret) #for CV folds and data splitting
library(ROCR) #for diagnostics and ROC plots/stats
library(gplots)
library(pROC) #same as ROCR
###Using 43 variables specified in the firstmodel data Appendix###
data.full<-data[,c("invasion.status", "es1", "es2", "es3","es4", "es5","es6","es7","es8","es9","es10","es11","es12","es13","es14","es15","es16","es17","es18","es19","es20","es21","es22","es23", "impg1", "impg2", "impn1", "impn2", "impn3", "impn4", "impn5", "impn6", "impa1", "impa2", "impa3", "impa4", "impp1", "impp2", "impp3", "impp4", "impp5", "impp6")]
###Converting DV into Factor with names for Caret Library###
data.full$invasion.status<-factor(data.full$invasion.status, levels=c(0,1), labels=c("MIN", "MAJ"))
set.seed(100) #the most metal seed for CV
#This method of data slicing-or CV will be used for logistic model
tc<-trainControl(method="repeatedcv", number=10, repeats = 10, summaryFunction=twoClassSummary, classProb=T)
#creates 10 CV folds for this data, summaryFunction provides a ROC summary stat in call to model#
###Logistic Regression Model Specification###
model.lr.1<-train(as.factor(invasion.status)~es1+es2+es3+es4+es5+es6+es7+es8+es9+es10+es11+es12+es13+es14+es15+es16+es17+es18+es19+es20+es21+es22+es23+impg1+impg2+impn1+impn2+impn3+impn4+impn5+impn6+impa1+impa2+impa3+impa4+impp1+impp2+impp3+impp4+impp5+impp6, metric="ROC", method="glm", family="binomial", trControl=tc, data=data.full)
summary(model.lr.1) #provides coefficients & traditional R model output
model.lr.1 #provides CV summary stats #keep in mind caret takes first class (here the first class is O) ##Generalized Linear Model 
#as reference class so sensitivity and specificity are backwards
confusionMatrix(model.lr.1, norm="average")
###Implementing RF (withCV) on entirety of data###
model.rf<-train(as.factor(invasion.status)~., metric="ROC", method="rf", importance=T, proximity=F, ntree=1000, trControl=tc, data=data.full)
#Variable importance measures are retained when importance=T
#Number of trees grown is indicated on ntree=1000
model.rf
confusionMatrix(model.rf, norm="average")
library(ROCR)
attach(data.full) #have to attach the data in order to get probabilities
###Gathering info for ROC Plots###
LR.1.pred<-predict(model.lr.1, data.full$invasion.status, type="prob")
RF.1.pred<-predict(model.rf, data.full$invasion.status, type="prob")
pred.LR<-prediction(LR.1.pred$MAJ, data$invasion.status)
perf.LR<-performance(pred.LR, "tpr", "fpr")
pred.RF.1<-prediction(RF.1.pred$MAJ, data$invasion.status)
perf.RF.1<-performance(pred.RF.1, "tpr", "fpr")
#plot(perf.LR, main="Logistic Regression and Random Forests")
#plot(perf.RF.1, add=T, lty=2)
#legend(0.32, 0.25, c("Logistic Regression 0.723", "Random Forest  0.885"), lty=c(1,2), bty="n", cex=.75)
par(mar = c(6.5, 6.5, 1.5, 0.5), mgp = c(5, 1, 0))
plot(perf.LR, main="Logistic Regression and Random Forests", lty=1, col=578)
plot(perf.RF.1, xlab = "False Positive Rate", ylab = "True Positive Rate", add=T, lty=3, col=84)
legend(0.50, 0.50, c("Logistic Regression 0.729", "Random Forest  0.877"), lty = c(1,3), col = c(578,84), bty="o", cex=.3)
dev.print(tiff, "secondmodel08-15-2017.tiff", res=600, height=5, width=7, units="in")
###Separation Plots###
library(separationplot)
##Have to transform DV back to 0,1 values for separation plots
data.full$invasion.status<-factor(data.full$invasion.status, levels=c(0,1), labels=c("MIN","MAJ"))
invasion.status<-as.vector(data.full$invasion.status) #transforming actual observations into vector
separationplot(LR.1.pred$MAJ, invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Logistic Regression", height=1.5, col0="white", col1="black")
separationplot(RF.1.pred$MAJ, invasion.status, type="line", line=T, lwd2=1, show.expected=T, heading="Random Forests", height=2.5, col0="white", col1="black")
